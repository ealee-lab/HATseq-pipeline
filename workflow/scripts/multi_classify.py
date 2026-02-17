import argparse
import shutil
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

def get_peaks_for_multiintvls(peak_list, intvls, donor):
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

def filter_multiintvls(peak_df, n):
    """Determine which intervals likely represent putative insertions. 
    These intervals will be re-classified based on information from
    multiple samples."""
    
    # Only re-classify putative insertions in >= n samples
    multi_peaks = peak_df[peak_df["num"] >= n]
    
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
        d=200, c=[5,6,7,9,10], 
        o=["collapse","collapse","collapse","collapse","max"])
    
    merged_df = merged_bt.to_dataframe(disable_auto_names=True, header=None)
    merged_df.columns = ["chrm","start","end","RPM","strand","classes","peaks","num"]
    
    merged_df.insert(3, "name", None)
    merged_df["name"] = donor + "-peak-" + merged_df["chrm"] + ":" + merged_df["start"].astype(str)

    return merged_df

def format_multiintvls(merged_df):
    """Re-format the merged intervals for subsequent steps."""
    format_df = merged_df.copy()
    
    format_df["RPM"] = format_df["RPM"].str.split(",")
    format_df["strand"] = format_df["strand"].str.split(",")
    format_df["classes"] = format_df["classes"].str.split(",")
    format_df["peaks"] = format_df["peaks"].str.split(",")

    format_df = format_df.explode(["RPM","strand","classes","peaks"])
    format_df = format_df.drop_duplicates(keep='first')

    format_df["RPM"] = format_df["RPM"].astype(float)
    format_df = format_df.groupby(["chrm","start","end","name","num"]).agg(
        RPM=('RPM', 'mean'),
        CoV=('RPM', 'std'),
        strand=('strand', lambda x: ",".join(set(x))),
        classes=('classes', lambda x: ",".join(x)),
        peaks=('peaks', lambda x: ",".join(set(x)))
    ).reset_index()
    format_df["CoV"] = format_df["CoV"] / format_df["RPM"]

    return format_df

def reclassify_peaks_by_label(format_df):
    """Re-classify peaks detected across multiple replicates."""
    reclass_df = format_df.copy()

    unk_counts = reclass_df["classes"].str.count("UNK")
    som_counts = reclass_df["classes"].str.count("SOM")
    private_counts = reclass_df["classes"].str.count("SOM_private")
    clonal_counts = reclass_df["classes"].str.count("SOM_clonal")

    # Classify as UNK if more reps are labeled UNK than SOM
    idxs = unk_counts > som_counts
    reclass_df.loc[idxs, "classification"] = "UNK"

    if private_counts.sum() == 0 & clonal_counts.sum() == 0:
        reclass_df.loc[~idxs, "classification"] = "SOM"
    else:
        # Classify as clonal if at least one rep is labeled clonal
        som_idxs = clonal_counts >= 1
        reclass_df.loc[som_idxs, "classification"] = "SOM_clonal"
        reclass_df.loc[~som_idxs, "classification"] = "SOM_private"

    # Peak is on different strands in different replicates
    reclass_df.loc[(reclass_df["strand"] == "+,-") |
                   (reclass_df["strand"] == "-,+"), "classification"] = "FP"

    return reclass_df

def reclassify_peaks_by_signal(format_df, n):
    """Re-classify putative insertions based on information from n tissues/cells."""
    reclass_df = format_df.copy()

    # Peak is in some samples
    reclass_df.loc[(reclass_df["num"] > 1) & 
                   (reclass_df["num"] < n), "classification"] = "SOM_clonal"
    
    # Peak is in all samples at similar RPMs
    reclass_df.loc[(reclass_df["num"] == n) &
                   (reclass_df["CoV"] < 0.5), "classification"] = "UNK"
    
    # Peak is in all samples with variable RPMs
    reclass_df.loc[(reclass_df["num"] == n) &
                   (reclass_df["CoV"] >= 0.5), "classification"] = "SOM_clonal"
    
    # Peak is on different strands in different samples
    reclass_df.loc[(reclass_df["strand"] == "+,-") |
                   (reclass_df["strand"] == "-,+"), "classification"] = "FP"
    
    return reclass_df

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
        all_df = all_df.merge(true_df, how="left", on="peak")

    all_df = all_df.drop_duplicates(keep='first')
    all_df.to_csv(peakfile, sep="\t", index=False, header=True)

def main(args):
    donor = args.donor
    filenames = args.filenames
    peakfile = args.peakfile
    multifile = args.multifile
    replicate = True if args.replicate == "rep" else False
    num_samples = len(filenames)
    min_samples = 2 if replicate else 3

    if num_samples < min_samples:
        # Concatenate files without comparison if too few samples
        dfs = []
        for file in filenames:
            file_df = pd.read_csv(file, sep="\t")
            dfs.append(file_df)
        df = pd.concat(dfs)
        df.to_csv(peakfile, sep="\t", index=False, header=True)
        
        # Create empty file
        with open(multifile, 'w') as f:
            pass

    else:
        dfs = []
        benchmark_dfs = []

        # Read in sample peak files
        for file in filenames:
            file_df = pd.read_csv(file, sep="\t")
            file_df = file_df[file_df["filter_reason"] != "nearby_peak"]
            file_df = file_df.sort_values(by=["chrm","start"])
        
            benchmark = True if "true_insertion_ID" in file_df.columns else False
            if benchmark:
                benchmark_dfs.append(file_df[["peak","true_insertion_ID"]])
                file_df = file_df.drop(columns="true_insertion_ID")

            dfs.append(file_df)
        
        # Save true insertion info if benchmarking
        if len(benchmark_dfs) != 0:
            benchmark_df = pd.concat(benchmark_dfs)

        bedtools = [BedTool.from_dataframe(df) for df in dfs]

        intvls = get_multiintvls(bedtools)
        peak_df = get_peaks_for_multiintvls(bedtools, intvls, donor)
        multi_df = filter_multiintvls(peak_df, min_samples)
        merge_df = merge_multiintvls(multi_df, donor)
        format_df = format_multiintvls(merge_df)

        if replicate:
            reclass_df = reclassify_peaks_by_label(format_df)
        else:
            reclass_df = reclassify_peaks_by_signal(format_df, num_samples)

        write_reclassified_peaks(reclass_df, dfs, peakfile, benchmark_df)
        write_multi_peaks(reclass_df, multifile)

    return

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
    parser.add_argument("-r", "--replicate", 
                        help="'rep' if comparing across replicates, else 'sample'",
                        required=True)
    args = parser.parse_args()
    main(args)
