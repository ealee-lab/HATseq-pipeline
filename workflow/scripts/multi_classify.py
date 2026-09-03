import os
import re
import argparse
import sys
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
from natsort import natsort_key

pd.options.mode.chained_assignment = None  # default='warn'

def read_sample_peaks(files):
    """Read in sample-specific peaks."""
    file_dfs = []
    analysis_dfs = []

    for file in files:
        if isinstance(file, (str, os.PathLike)):
            file_df = pd.read_csv(file, sep="\t")
            file_df.columns = file_df.columns.str.strip("#")
        elif isinstance(file, pd.DataFrame):
            file_df = file
            file_df.columns = file_df.columns.str.strip("#")
        else:
            raise ValueError("Invalid data type")
        
        analysis_df = file_df.iloc[:, 0:8]

        # Remove known peaks
        analysis_df = analysis_df[
            ~analysis_df["classification"].str.contains(r'KR|KNR|Non-spec')]

        # Remove artifacts and peaks in error-prone regions
        analysis_df = analysis_df[
            ~analysis_df["filter"].str.contains(r'nearby|error')]

        # Sort peaks for downstream analysis
        analysis_df = analysis_df.sort_values(by=["chrm","start"])

        file_dfs.append(file_df)
        analysis_dfs.append(analysis_df)
    
    return analysis_dfs, file_dfs

def get_multiintvls(peak_BTs):
    """Find all covered intervals, including which samples intersect each interval."""
    x = BedTool()
    fns = [bt.fn for bt in peak_BTs]
    intvls_BT = x.multi_intersect(i=fns)
    intvls_BT = intvls_BT.cut([0,1,2,4])
    print(intvls_BT)
    intvls_BT = intvls_BT.sort().merge(d=50, c=4, o="distinct")
    return intvls_BT

def get_peaks_for_multiintvls(peak_BTs, intvls_BT, donor):
    """Find the corresponding sample-specific peak for each interval."""
    sample_dfs = []

    # Intersect sample-specific peaks with intervals
    for sample_BT in peak_BTs:
        sample_intvls = sample_BT.intersect(intvls_BT, wa=True, wb=True)
        sample_df = sample_intvls.to_dataframe(disable_auto_names=True, header=None)
        
        sample_df.columns = [
            "sample_chrm","sample_start","sample_end","peak",
            "RPM","strand","classification","filter",
            "chrm","start","end","file_nums"]
        
        sample_df["num_samples"] = sample_df["file_nums"].apply(lambda x: len(set(x.split(","))))
        sample_dfs.append(sample_df)
    
    peak_df = pd.concat(sample_dfs)

    # Add common name for peaks belonging to same interval
    peak_df["name"] = donor + "-peak-" + peak_df["chrm"] + ":" + peak_df["start"].astype(str)
    peak_df = peak_df[["chrm","start","end","name","RPM","strand",
                       "classification","filter","peak","num_samples"]]

    # peak_df = peak_df.sort_values(by="name")
    return peak_df

def filter_multiintvls(peak_df):
    """Filter for intervals that likely represent putative insertions."""
    # Only re-classify putative insertions in >= 2 samples
    multi_peaks = peak_df[peak_df["num_samples"] >= 2]
    
    # Interval must pass filters in >= 1 sample
    put_peaks = multi_peaks[multi_peaks["filter"] == "PASS"]["name"].unique()
    multi_df = multi_peaks[multi_peaks["name"].isin(put_peaks)]
    return multi_df

def format_multiintvls(multi_df):
    """Re-format the intervals for subsequent steps."""
    format_df = multi_df.drop_duplicates(keep='first')

    format_df = format_df.groupby(["chrm","start","end","name","num_samples"]).agg(
        RPM=('RPM', lambda x: ",".join(x.astype(str))),
        strand=('strand', lambda x: ",".join(set(x))),
        classes=('classification', lambda x: ",".join(x)),
        peaks=('peak', lambda x: ",".join(set(x)))
    ).reset_index()
    format_df["filter"] = "."

    return format_df

def reclassify_peaks_by_label(format_df):
    """Re-classify peaks detected across multiple replicates."""
    reclass_df = format_df.copy()

    unk_counts = reclass_df["classes"].str.count("UNK")
    som_counts = reclass_df["classes"].str.count("SOM")
    private_counts = reclass_df["classes"].str.count("SOM_private")
    clonal_counts = reclass_df["classes"].str.count("SOM_clonal")

    # Call UNK if more reps are labeled UNK than SOM
    unk_idxs = unk_counts > som_counts
    reclass_df.loc[unk_idxs, "classification"] = "UNK"
    reclass_df.loc[~unk_idxs, "classification"] = "SOM"

    # Where there is a distinction between clonal and private peaks...
    if private_counts.sum() > 0 or clonal_counts.sum() > 0:
        clonal_idxs = clonal_counts >= 1 # call clonal if >= 1 rep is labeled clonal
        reclass_df.loc[~unk_idxs & clonal_idxs, "classification"] = "SOM_clonal"
        reclass_df.loc[~unk_idxs & ~clonal_idxs, "classification"] = "SOM_private"
    
    # Update filter for all re-classified peaks
    reclass_df["filter"] = "PASS (multi-sample)"

    # Call FP if peak is on different strands in different replicates
    reclass_df.loc[(reclass_df["strand"] == "+,-") |
                   (reclass_df["strand"] == "-,+"), "classification"] = "FP"
    reclass_df.loc[reclass_df["classification"] == "FP", "filter"] = "strand (multi-sample)"

    return reclass_df

def reclassify_peaks_by_signal(format_df, total_samples):
    """Re-classify peaks detected across multiple tissues/cells."""
    reclass_df = format_df.copy()

    # UNK threshold is 2/3 total number of samples
    thresh = np.ceil((2 * total_samples) / 3)

    # TODO: Combine RPM for replicates before this step so that all_RPM cannot be greater than threshold
    reclass_df["all_RPM"] = reclass_df["RPM"].str.split(",").apply(
        lambda x: sum([float(n) >= 25 for n in x]))
    
    # Peak is in some samples
    reclass_df.loc[(reclass_df["num_samples"] > 1) & 
                   (reclass_df["num_samples"] < thresh), "classification"] = "SOM_clonal"
    
    # Peak is in most samples with high RPM
    reclass_df.loc[(reclass_df["num_samples"] >= thresh) &
                   (reclass_df["all_RPM"] >= thresh), "classification"] = "UNK"
    
    # Peak is in most samples with low RPM
    reclass_df.loc[(reclass_df["num_samples"] >= thresh) &
                   (reclass_df["all_RPM"] < thresh), "classification"] = "SOM_clonal"
    
    # Update filter for all re-classified peaks
    reclass_df["filter"] = "PASS (multi-sample)"

    # Peak is on different strands in different samples
    reclass_df.loc[(reclass_df["strand"] == "+,-") |
                   (reclass_df["strand"] == "-,+"), "classification"] = "FP"
    reclass_df.loc[reclass_df["classification"] == "FP", "filter"] = "strand (multi-sample)"
    
    return reclass_df

def reclassify_KNR_peaks(format_df, total_samples):
    """Re-classify peaks detected across multiple samples."""
    reclass_df = format_df.copy()

    thresh = np.ceil(total_samples / 3)
    
    # Pass if KNR is in at least a third of samples
    reclass_df.loc[(reclass_df["num_samples"] >= thresh), "classification"] = "KNR"
    reclass_df["filter"] = "PASS (multi-sample)"

    # Set to lowConf if KNR is in less than a third of samples
    reclass_df.loc[reclass_df["num_samples"] < thresh, "filter"] = "lowConf (multi-sample)"

    return reclass_df

def get_multi_peaks(reclass_df):
    """Get all re-classified multi-sample peaks for the given donor."""
    out_df = reclass_df[["chrm","start","end","name","RPM","strand",
                         "classification","peaks","num_samples"]]
    out_df = out_df.sort_values(by="name", key=natsort_key)
    return out_df

def get_knr_peaks(reclass_df):
    """Get final KNR list for the given donor."""
    out_df = reclass_df[["chrm","start","end","name","RPM","strand"]]
    # out_df["RPM"] = out_df["RPM"].str.split(",").apply(lambda x: np.median([float(n) for n in x]))
    out_df = out_df.sort_values(by="name", key=natsort_key)
    return out_df

def get_reclassified_peaks(reclass_df, peaks_list):
    """Get all re-classified sample-specific peaks for the given donor."""
    col_df = reclass_df[["classification","filter","peaks"]]
    col_df["peaks"] = col_df["peaks"].str.split(",")
    col_df = col_df.explode("peaks").reset_index(drop=True)
    col_df.columns = ["donor_class","donor_filter","peak"]

    out_dfs = []
    for sample_peaks in peaks_list:
        out_df = sample_peaks.merge(col_df, how="left", on="peak")
        out_df.loc[~out_df["donor_class"].isna(), "classification"] = out_df["donor_class"]
        out_df.loc[~out_df["donor_class"].isna(), "filter"] = out_df["donor_filter"]
        out_df.drop(columns=["donor_class","donor_filter"], inplace=True)
        out_dfs.append(out_df)

    all_df = pd.concat(out_dfs).sort_values(by=["chrm","start"], key=natsort_key)

    all_df = all_df.drop_duplicates(keep='first')
    return all_df

def process_peaks(analysis_dfs, donor, comparison):
    """Given a list of sample-specific peaks, find and format intervals where
    the peaks overlap across multiple samples."""
    peak_BTs = [BedTool.from_dataframe(df) for df in analysis_dfs]

    intvls_BT = get_multiintvls(peak_BTs)
    peak_df = get_peaks_for_multiintvls(peak_BTs, intvls_BT, donor)
    multi_df = filter_multiintvls(peak_df)
    format_df = format_multiintvls(multi_df)

    return format_df

def run_comparison(donor, peaks, comparison, total_samples):
    """Run the multi-sample comparison for a given donor."""
    # Compare SOM and UNK peaks across samples
    analysis_dfs, file_dfs = read_sample_peaks(peaks)
    format_som_df = process_peaks(analysis_dfs, donor, comparison)

    if comparison == "rep":
        reclass_som_df = reclassify_peaks_by_label(format_som_df)
    else:
        reclass_som_df = reclassify_peaks_by_signal(format_som_df, total_samples)

    # Compare KNR peaks across samples
    knr_dfs = [df[df["classification"] == "KNR"] for df in file_dfs]
    
    format_knr_df = process_peaks(knr_dfs, donor, comparison)
    reclass_knr_df = reclassify_KNR_peaks(format_knr_df, total_samples)

    # Output dataframes
    peakfile_df = get_reclassified_peaks(pd.concat([reclass_som_df, reclass_knr_df]), file_dfs)
    multifile_df = get_multi_peaks(reclass_som_df)
    knownfile_df = get_knr_peaks(reclass_knr_df)

    return peakfile_df, multifile_df, knownfile_df

def main(donor, filenames, peakfile, multifile, knownfile, comparison):
    # Cross-replicate OR cross-tissue/cell comparison
    if comparison != "both":
        total_samples = len(filenames)
        
        out_df1, out_df2, out_df3 = run_comparison(
            donor, filenames, comparison, total_samples)
        
        out_df1.to_csv(peakfile, sep="\t", index=False, header=True)
        out_df2.to_csv(multifile, sep="\t", index=False, header=True)
        out_df3.to_csv(knownfile, sep="\t", index=False, header=True)

    # Cross-replicate AND cross-tissue comparison
    else:
        tissue_files = {}

        for file in filenames:
            tissue = re.search(rf'.*/{donor}-(.*)_[A-Z][0-9]+_classified_peaks\.bed', file).group(1)

            if tissue not in tissue_files:
                tissue_files[tissue] = [file]
            else:
                tissue_files[tissue].append(file)
        
        peakfiles = []
        multifiles = []
        knownfiles = []

        for tissue in tissue_files:
            num_reps = len(tissue_files[tissue])

            if num_reps < 2:
                tissue_df = pd.read_csv(tissue_files[tissue][0], sep="\t")
                tissue_df.columns = tissue_df.columns.str.strip("#")
                peakfiles.append(tissue_df)
            else:
                tissue_df1, tissue_df2, tissue_df3 = run_comparison(
                    donor, tissue_files[tissue], "rep", num_reps)
                peakfiles.append(tissue_df1)

        num_tissues = len(tissue_files)
        if num_tissues < 2:
            raise ValueError("At least two tissues are required to use 'both' comparison.")
        else:
            both_df1, both_df2, both_df3 = run_comparison(donor, peakfiles, "tissue", num_tissues)
            both_df1.to_csv(peakfile, sep="\t", index=False, header=True)
            both_df2.to_csv(multifile, sep="\t", index=False, header=True)
            both_df3.to_csv(knownfile, sep="\t", index=False, header=True)

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="A script to classify peaks present in multiple samples for a given donor.")
    
    parser.add_argument("-d", "--donor", 
                        help="Donor name", 
                        required=True)
    parser.add_argument("-f", "--filenames", 
                        help="List of file names of all sample-specific peaks for donor",
                        nargs="+", 
                        required=True)
    parser.add_argument("-p", "--peakfile", 
                        help="Path to output file with all (re-)classified peaks for donor",
                        required=True)
    parser.add_argument("-m", "--multifile", 
                        help="Path to output file with multi-sample peaks for donor",
                        required=True)
    parser.add_argument("-k", "--knownfile", 
                        help="Path to output file with final KNR peaks for donor",
                        required=True)
    parser.add_argument("-c", "--comparison", 
                        help="Mode of comparison (tissue, cell, rep, or both)",
                        required=True)
    parser.add_argument("-l", "--log", 
                        help="Path to log file",
                        required=True)
    args = parser.parse_args()

    sys.stderr = open(args.log, "w", buffering=1)

    if len(args.filenames) < 2:
        raise ValueError("At least two sample-specific peak files are required for comparison.")
    else:
        main(args.donor, args.filenames, args.peakfile, 
             args.multifile, args.knownfile, args.comparison)
    