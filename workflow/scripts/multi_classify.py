import argparse
import pandas as pd
import pybedtools
from pybedtools import BedTool

def get_multiintvls(names):
    """Report all intervals along with the number of samples that 
    intersect each interval."""
    x = BedTool()
    intvls = x.multi_intersect(i=names)
    intvls = intvls.merge(c=4, o='max')
    return intvls

def get_corresponding_sample_peaks(peaks, intvls, donor):
    """Find the corresponding tissue-specific peak for each interval."""
    dfs = []

    for sample in peaks:
        sample_intvls = sample.intersect(intvls, wa=True, wb=True)
        df = sample_intvls.to_dataframe(disable_auto_names=True, header=None)
        df.columns = ["s_chrm","s_start","s_end","peak",
                      "RPM","strand","classification",
                      "chrm","start","end","num"]
        dfs.append(df)
    peak_df = pd.concat(dfs)

    peak_df["name"] = donor + "-peak-" + peak_df["chrm"] + ":" + peak_df["start"].astype(str)
    peak_df = peak_df.sort_values(by="name")

    return peak_df

def get_peaks_for_reclassification(peak_df):
    """Determine which peaks likely represent putative insertions, and thus, 
    which can be re-classified based on multi-tissue information."""
    # Only re-classify putative insertions in >= 3 tissues
    multi_peaks = peak_df[peak_df["num"] >= 3]
    put_peaks = multi_peaks[multi_peaks["classification"].isin(["UNK","SOM_clonal","SOM_private"])]["name"].unique()
    known_peaks = multi_peaks[multi_peaks["classification"].isin(["KR","KNR","Off-target"])]["name"].unique()

    # Peak must be UNK/SOM in >= 1 tissue and not KR/KNR/Off-target in any tissue
    reclass_df = multi_peaks[(multi_peaks["name"].isin(put_peaks)) & (~multi_peaks["name"].isin(known_peaks))]
    peak_map = dict(zip(reclass_df["peak"], reclass_df["name"]))

    group_df = reclass_df.groupby(["chrm","start","end","name","num"])
    group_df = group_df.agg(
        RPM=pd.NamedAgg("RPM", lambda x: ",".join(str(x))),
        CV=pd.NamedAgg("RPM", lambda x: x.std() / x.mean()),
        strand=pd.NamedAgg("strand", lambda x: ",".join(set(x))),
        peaks=pd.NamedAgg("peak", lambda x: ",".join(x)))
    
    return group_df.reset_index(), peak_map

def reclassify_peaks(reclass_df, num):
    """Re-classify putative insertions based on information from multiple tissues."""
    reclass_df["classification"] = "-NA-"
    
    reclass_df.loc[(reclass_df["num"] > 1) & 
                   (reclass_df["num"] < num), "classification"] = "SOM_clonal"
    
    reclass_df.loc[(reclass_df["num"] == num) &
                   (reclass_df["CV"] < 0.5), "classification"] = "UNK"
    reclass_df.loc[(reclass_df["num"] == num) &
                   (reclass_df["CV"] >= 0.5), "classification"] = "SOM_clonal"
    
    reclass_df.loc[(reclass_df["strand"] == "+,-") |
                   (reclass_df["strand"] == "-,+"), "classification"] = "FP"
    
    return

def format_and_write_peaks(peaks, reclass_df, peak_map, peakfile, multifile):
    """Concatenate and write all relevant peaks for the given donor."""
    all_df = pd.concat(peaks)
    put_df = all_df[all_df["classification"].isin(["FP","UNK","SOM_clonal","SOM_private"])]

    put_df["donor_peak"] = put_df["peak"].map(peak_map)
    put_df.to_csv(peakfile, sep="\t", index=False, header=True)

    multi_df = reclass_df[["chrm","start","end","name","num","strand","classification","peaks"]]
    multi_df.to_csv(multifile, sep="\t", index=False, header=True)

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to classify peaks occurring in multiple tissues.")
    parser.add_argument("-f", "--filenames", 
                        help="List of file names for all sample peaks for given donor", 
                        nargs="+", 
                        required=True)
    parser.add_argument("-n", "--num", 
                        help="Number of tissues for given donor", 
                        type=int,
                        required=True)
    parser.add_argument("-d", "--donor", 
                        help="Donor name", 
                        required=True)
    parser.add_argument("-p", "--peakfile", 
                        help="Path to output file with all putative peaks for donor",
                        required=True)
    parser.add_argument("-m", "--multifile", 
                        help="Path to output file with multi-tissue peaks for donor",
                        required=True)
    args = parser.parse_args()

    # Read in sample peak files
    dfs = []
    for file in args.filenames:
        df = pd.read_csv(file, sep="\t", 
                         usecols=[0,1,2,3,4,5,6], 
                         names=["chrm","start","end","peak","RPM","strand","classification"])
        dfs.append(df)
    bedtools = [BedTool.from_dataframe(df) for df in dfs]

    intvls = get_multiintvls(args.filenames)
    peak_df = get_corresponding_sample_peaks(bedtools, intvls, args.donor)
    reclass_df, peak_map = get_peaks_for_reclassification(peak_df)
    reclassify_peaks(reclass_df, args.num)
    format_and_write_peaks(dfs, reclass_df, peak_map, args.peakfile, args.multifile)
