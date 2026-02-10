import argparse
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
from natsort import natsort_key

def get_multiintvls(bedtools):
    """Report all covered intervals along with the number of samples that
    intersect each interval."""
    x = BedTool()
    fns = [bt.fn for bt in bedtools]
    intvls = x.multi_intersect(i=fns)
    intvls = intvls.cut([0,1,2,3])
    return intvls

def get_sample_peaks_for_multiintvls(peak_list, intvls, donor):
    """Find the corresponding sample-specific peak for each interval."""
    sample_dfs = []

    for sample_peaks in peak_list:
        sample_intvls = sample_peaks.intersect(intvls, wa=True, wb=True)
        sample_df = sample_intvls.to_dataframe(disable_auto_names=True, header=None)

        sample_df.columns = [
            "chrm_s","start_s","end_s","peak",
            "RPM","strand","classification","filter",
            "chrm","start","end","num"]
        sample_dfs.append(sample_df)
    
    peak_df = pd.concat(sample_dfs)

    peak_df["name"] = donor + "-peak-" + peak_df["chrm"] + ":" + peak_df["start"].astype(str)
    peak_df = peak_df[["chrm","start","end","name","RPM","strand",
                       "classification","filter","peak","num"]]

    peak_df = peak_df.sort_values(by="name")

    return peak_df

def filter_multiintvls(peak_df):
    """Determine which intervals likely represent putative insertions. 
    These intervals will be re-classified based on information from
    multiple samples."""
    
    # Only re-classify putative insertions in >= 3 samples
    multi_peaks = peak_df[peak_df["num"] >= 3]
    
    put_peaks = multi_peaks[(multi_peaks["classification"] == "UNK") | 
                            (multi_peaks["classification"].str.contains("SOM"))]["name"].unique()
    known_peaks = multi_peaks[multi_peaks["classification"].isin(["KR","KNR","Off-target"])]["name"].unique()
    filter_peaks = multi_peaks[~(multi_peaks["filter"].isin(["-NA-","RPM","polyApercent","templates"]))]["name"].unique()

    # Interval must be UNK/SOM in >= 1 sample... 
    # AND not KR/KNR/Off-target in any sample...
    # AND not labeled as artifact in any sample
    multi_df = multi_peaks[multi_peaks["name"].isin(put_peaks)]
    multi_df = multi_df[~multi_df["name"].isin(known_peaks)]
    multi_df = multi_df[~multi_df["name"].isin(filter_peaks)]

    return multi_df

def merge_multiintvls(multi_df, donor):
    """Merge adjacent intervals."""
    sort_df = multi_df.sort_values(by=["chrm","start"])
    merged_bt = BedTool.from_dataframe(sort_df).merge(
        d=50, c=[5,6,9,10], o=["collapse","distinct","distinct","max"])
    
    merged_df = merged_bt.to_dataframe(disable_auto_names=True, header=None)
    merged_df.columns = ["chrm","start","end","RPM","strand","peaks","num"]
    
    merged_df.insert(3, "name", None)
    merged_df["name"] = donor + "-peak-" + merged_df["chrm"] + ":" + merged_df["start"].astype(str)

    return merged_df

def reclassify_peaks(merged_df, n):
    """Re-classify putative insertions based on information from n samples."""
    # Find variation in signal for each peak
    RPM_std = merged_df["RPM"].apply(lambda x: np.std([float(n) for n in x.split(",")]))
    RPM_mean = merged_df["RPM"].apply(lambda x: np.mean([float(n) for n in x.split(",")]))
    
    reclass_df = merged_df.copy()
    reclass_df["RPM"] = RPM_mean
    reclass_df["CV"] = RPM_std / RPM_mean

    ## Re-classify peaks ##
    reclass_df["classification"] = "-NA-"

    # Peak is in some samples
    reclass_df.loc[(reclass_df["num"] > 1) & 
                   (reclass_df["num"] < n), "classification"] = "SOM_clonal"
    
    # Peak is in all samples at similar RPMs
    reclass_df.loc[(reclass_df["num"] == n) &
                   (reclass_df["CV"] < 0.5), "classification"] = "UNK"
    
    # Peak is in all samples with variable RPMs
    reclass_df.loc[(reclass_df["num"] == n) &
                   (reclass_df["CV"] >= 0.5), "classification"] = "SOM_clonal"
    
    # Peak is on different strands in different samples
    reclass_df.loc[(reclass_df["strand"] == "+,-") |
                   (reclass_df["strand"] == "-,+"), "classification"] = "FP"
    
    return reclass_df

def write_candidate_peaks(reclass_df, peak_list, peakfile):
    """Concatenate and write all peaks for the given donor to 
    one file."""
    all_df = pd.concat(peak_list)
    # cand_df = all_df[(all_df["classification"] == "UNK") |
    #                  (all_df["classification"].str.contains("SOM")) |
    #                  (all_df["classification"] == "FP")]

    donor_df = reclass_df[["name","peaks","classification"]]
    
    donor_df.loc[:, "peaks"] = donor_df["peaks"].str.split(",")
    donor_df = donor_df.explode("peaks").reset_index(drop=True)
    peak_map = dict(zip(donor_df["peaks"], donor_df["name"]))
    all_df.loc[:, "donor_peak"] = all_df["peak"].map(peak_map)

    donor_df = donor_df.drop("peaks", axis=1)
    donor_df.columns = ["donor_peak","donor_classification"]
    out_df = all_df.merge(donor_df, how="left", on="donor_peak")
    
    out_df = out_df.drop_duplicates(keep="first")
    out_df = out_df.sort_values(by=["chrm","start"], key=natsort_key)
    out_df.to_csv(peakfile, sep="\t", index=False, header=True)

def write_multi_peaks(reclass_df, multifile):
    """Write all re-classified multi-sample peaks for the given donor."""
    out_df = reclass_df[["chrm","start","end","name","RPM","strand","classification","peaks","num"]]
    out_df = out_df.sort_values(by="name", key=natsort_key)
    out_df.to_csv(multifile, sep="\t", index=False, header=True)

def write_reclassified_peaks(reclass_df, peak_list, peakfile, true_df=None):
    col_df = reclass_df[["classification","peaks"]]
    col_df["peaks"] = col_df["peaks"].str.split(",")
    col_df = col_df.explode("peaks").reset_index(drop=True)
    col_df.columns = ["donor_class","peak"]

    out_dfs = []
    for sample_peaks in peak_list:
        out_df = sample_peaks.merge(col_df, how="left", on="peak")
        out_df.loc[~out_df["donor_class"].isna(), "classification"] = out_df["donor_class"]
        out_df.loc[~out_df["donor_class"].isna(), "filter_reason"] = "multi-sample"
        out_df.drop(columns="donor_class", inplace=True)
        out_dfs.append(out_df)

    all_df = pd.concat(out_dfs)

    if isinstance(true_df, pd.DataFrame):
        all_df = all_df.merge(benchmark_df, how="left", on="peak")

    all_df.to_csv(peakfile, sep="\t", index=False, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to classify peaks occurring in multiple samples.")
    
    parser.add_argument("-d", "--donor", 
                        help="Donor name", 
                        required=True)
    parser.add_argument("-f", "--filenames", 
                        help="List of file names for all sample peaks for given donor", 
                        nargs="+", 
                        required=True)
    parser.add_argument("-p", "--peakfile", 
                        help="Path to output file with all putative peaks for donor",
                        required=True)
    parser.add_argument("-m", "--multifile", 
                        help="Path to output file with multi-sample peaks for donor",
                        required=True)
    args = parser.parse_args()

    num_samples = len(args.filenames)

    peak_dfs = []
    benchmark_dfs = []

    # Read in sample peak files
    for file in args.filenames:
        file_df = pd.read_csv(file, sep="\t")
        file_df = file_df[file_df["filter_reason"] != "nearby_peak"]
        file_df = file_df.sort_values(by=["chrm","start"])

        benchmark = True if "true_insertion_ID" in file_df.columns else False
        if benchmark:
            benchmark_dfs.append(file_df[["peak","true_insertion_ID"]])
            file_df = file_df.drop(columns="true_insertion_ID")

        peak_dfs.append(file_df)
    
    # Save true insertion info if benchmarking
    if len(benchmark_dfs) != 0:
        benchmark_df = pd.concat(benchmark_dfs)

    peak_bts = [BedTool.from_dataframe(df) for df in peak_dfs]

    intervals = get_multiintvls(peak_bts)
    interval_peaks = get_sample_peaks_for_multiintvls(peak_bts, intervals, args.donor)
    filtered_intervals = filter_multiintvls(interval_peaks)
    merged_intervals = merge_multiintvls(filtered_intervals, args.donor)
    reclassed_peaks = reclassify_peaks(merged_intervals, num_samples)
    write_reclassified_peaks(reclassed_peaks, peak_dfs, args.peakfile, benchmark_df)
    write_multi_peaks(reclassed_peaks, args.multifile)
