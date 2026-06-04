import os
import re
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

def get_peaks_for_multiintvls(peaks_list, intvls, donor):
    """Find the corresponding sample-specific peak for each interval."""
    sample_dfs = []

    for sample_peaks in peaks_list:
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
    """Filter for intervals that likely represent putative insertions."""
    # Only re-classify putative insertions in >= n samples
    multi_peaks = peak_df[peak_df["num"] >= n]
    
    put_peaks = multi_peaks[(multi_peaks["classification"] == "UNK") | 
                            (multi_peaks["classification"].str.contains("SOM"))]["name"].unique()

    # Interval must be UNK/SOM in >= 1 sample
    multi_df = multi_peaks[multi_peaks["name"].isin(put_peaks)]
    
    return multi_df

def merge_multiintvls(multi_df, donor):
    """Merge adjacent intervals."""
    sort_df = multi_df.sort_values(by=["chrm","start"])
    merged_bt = BedTool.from_dataframe(sort_df).merge(
        d=150, c=[5,6,7,9,10], 
        o=["collapse","collapse","collapse","collapse","max"])
    
    merged_df = merged_bt.to_dataframe(disable_auto_names=True, header=None)
    merged_df.columns = ["chrm","start","end","RPM","strand","classes","peaks","num"]
    
    merged_df.insert(3, "name", None)
    merged_df["name"] = donor + "-peak-" + merged_df["chrm"] + ":" + merged_df["start"].astype(str)

    return merged_df

def extract_num(peaks, comparison):
    """Re-count the number of samples corresponding to each interval.
    Merging may join multiple sub-intervals, necessitating a re-count."""
    if comparison == "rep":
        # Extract replicate substring
        samples = peaks.str.extract(r'.*_([A-Z][0-9]+)-(?:plus|minus)-peak-[0-9]+')
    else:
        # Extract tissue substring
        # Case 1: Peak has tissue AND replicate info
        # Case 2: No match, peak only has tissue info
        samples = peaks.str.extract(r'(.*)_[A-Z][0-9]+-(?:plus|minus)-peak-[0-9]+').fillna(
            peaks.str.extract(r'(.*)-(?:plus|minus)-peak-[0-9]+')
        )

    return samples[0].nunique()

def format_multiintvls(merged_df, comparison):
    """Re-format the merged intervals for subsequent steps."""
    format_df = merged_df.copy()
    
    format_df["RPM"] = format_df["RPM"].str.split(",")
    format_df["strand"] = format_df["strand"].str.split(",")
    format_df["classes"] = format_df["classes"].str.split(",")
    format_df["peaks"] = format_df["peaks"].str.split(",")

    format_df = format_df.explode(["RPM","strand","classes","peaks"])
    format_df = format_df.drop_duplicates(keep='first')

    format_df["RPM"] = format_df["RPM"].astype(float)
    format_df = format_df.groupby(["chrm","start","end","name"]).agg(
        RPM=('RPM', 'mean'),
        CoV=('RPM', 'std'),
        strand=('strand', lambda x: ",".join(set(x))),
        classes=('classes', lambda x: ",".join(x)),
        peaks=('peaks', lambda x: ",".join(set(x))),
        num=('peaks', lambda x: extract_num(x, comparison))
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
    unk_idxs = unk_counts > som_counts
    reclass_df.loc[unk_idxs, "classification"] = "UNK"
    reclass_df.loc[~unk_idxs, "classification"] = "SOM"

    # Where there is a distinction between clonal and private peaks...
    if private_counts.sum() > 0 or clonal_counts.sum() > 0:
        clonal_idxs = clonal_counts >= 1 # call clonal if >= 1 rep is labeled clonal
        reclass_df.loc[~unk_idxs & clonal_idxs, "classification"] = "SOM_clonal"
        reclass_df.loc[~unk_idxs & ~clonal_idxs, "classification"] = "SOM_private"

    # Peak is on different strands in different replicates
    reclass_df.loc[(reclass_df["strand"] == "+,-") |
                   (reclass_df["strand"] == "-,+"), "classification"] = "FP"

    return reclass_df

def reclassify_peaks_by_signal(format_df, n):
    """Re-classify peaks based on information from n tissues/cells."""
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

def get_multi_peaks(reclass_df):
    """Get all re-classified multi-sample peaks for the given donor."""
    out_df = reclass_df[["chrm","start","end","name","RPM","strand","classification","peaks","num"]]
    out_df = out_df.sort_values(by="name", key=natsort_key)
    return out_df

def get_reclassified_peaks(reclass_df, peaks_list):
    col_df = reclass_df[["classification","peaks"]]
    col_df["peaks"] = col_df["peaks"].str.split(",")
    col_df = col_df.explode("peaks").reset_index(drop=True)
    col_df.columns = ["donor_class","peak"]

    out_dfs = []
    for sample_peaks in peaks_list:
        out_df = sample_peaks.merge(col_df, how="left", on="peak")
        out_df.loc[~out_df["donor_class"].isna(), "classification"] = out_df["donor_class"]
        out_df.loc[~out_df["donor_class"].isna(), "filter_reason"] = "multi-sample"
        out_df.drop(columns="donor_class", inplace=True)
        out_dfs.append(out_df)

    all_df = pd.concat(out_dfs).sort_values(by=["chrm","start"], key=natsort_key)

    # if len(true_df) > 0:
    #     all_df = all_df.merge(true_df, how="left", on="peak")

    all_df = all_df.drop_duplicates(keep='first')
    return all_df

def get_cat_peaks(filenames):
    """Concatenate peaks given file names."""
    file_dfs = []
    for file in filenames:
        file_df = pd.read_csv(file, sep="\t")
        file_dfs.append(file_df)
    df = pd.concat(file_dfs)

    return df, pd.DataFrame()

def read_sample_peaks(files):
    file_dfs = []
    analysis_dfs = []

    for file in files:
        if isinstance(file, (str, os.PathLike)):
            file_df = pd.read_csv(file, sep="\t")
        elif isinstance(file, pd.DataFrame):
            file_df = file
        else:
            raise ValueError("Invalid data type")
        
        analysis_df = file_df.iloc[:, 0:8]

        # Remove known peaks
        analysis_df = analysis_df[~analysis_df["classification"].str.contains(r'KR|KNR|Non-specific')]

        # Remove artifacts and peaks in error-prone regions
        analysis_df = analysis_df[~analysis_df["filter_reason"].str.contains(r'nearby|error-prone')]

        analysis_df = analysis_df.sort_values(by=["chrm","start"])

        file_dfs.append(file_df)
        analysis_dfs.append(analysis_df)
    
    return analysis_dfs, file_dfs

def run_comparison(peaks, donor, comparison, num_samples, min_samples):
    if num_samples < min_samples:
        if isinstance(peaks[0], (str, os.PathLike)):
            # Concatenate files without comparison if too few samples
            peakfile_df, multifile_df = get_cat_peaks(peaks)
        else:
            raise ValueError("Must use file names")
    else:
        analysis_dfs, file_dfs = read_sample_peaks(peaks)

        bedtools = [BedTool.from_dataframe(df) for df in analysis_dfs]
        intvls = get_multiintvls(bedtools)
        peak_df = get_peaks_for_multiintvls(bedtools, intvls, donor)
        multi_df = filter_multiintvls(peak_df, min_samples)
        merged_df = merge_multiintvls(multi_df, donor)
        format_df = format_multiintvls(merged_df, comparison)

        if comparison == "rep":
            reclass_df = reclassify_peaks_by_label(format_df)
        else:
            reclass_df = reclassify_peaks_by_signal(format_df, num_samples)

        peakfile_df = get_reclassified_peaks(reclass_df, file_dfs)
        multifile_df = get_multi_peaks(reclass_df)

    return peakfile_df, multifile_df

def main(args):
    donor = args.donor
    filenames = args.filenames
    peakfile = args.peakfile
    multifile = args.multifile
    comparison = args.comparison

    if comparison != "both":
        num_samples = len(filenames)
        min_samples = 2 if comparison == "rep" else 3

        out_df1, out_df2 = run_comparison(filenames, donor, comparison, 
                                          num_samples, min_samples)
        out_df1.to_csv(peakfile, sep="\t", index=False, header=True)
        out_df2.to_csv(multifile, sep="\t", index=False, header=True)

    else:
        tissue_files = {}

        for file in filenames:
            tissue = re.search(f'.*/{donor}-(.*)_[A-Z][0-9]+_.*', file).group(1)

            if tissue not in tissue_files:
                tissue_files[tissue] = [file]
            else:
                tissue_files[tissue].append(file)
        
        peakfiles = []
        multifiles = []

        for tissue in tissue_files:
            num_reps = len(tissue_files[tissue])
            tissue_df1, tissue_df2 = run_comparison(tissue_files[tissue], donor, "rep", 
                                                    num_reps, 2)
            peakfiles.append(tissue_df1)
            multifiles.append(tissue_df2)

        num_tissues = len(tissue_files)
        if num_tissues < 3:
            out_df1 = pd.concat(peakfiles)
            out_df2 = pd.concat(multifiles)
            
            out_df1.to_csv(peakfile, sep="\t", index=False, header=True)
            out_df2.to_csv(multifile, sep="\t", index=False, header=True)
        else:
            both_df1, both_df2 = run_comparison(peakfiles, donor, "tissue", 
                                                num_tissues, 3)
            both_df1.to_csv(peakfile, sep="\t", index=False, header=True)
            both_df2.to_csv(multifile, sep="\t", index=False, header=True)

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
    parser.add_argument("-c", "--comparison", 
                        help="How to compare samples (tissue, cell, rep, or both)",
                        required=True)
    args = parser.parse_args()
    main(args)
