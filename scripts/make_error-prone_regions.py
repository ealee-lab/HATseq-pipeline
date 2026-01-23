import os
import glob
import argparse
import pandas as pd
import pybedtools
from pybedtools import BedTool

# def concat_bed_by_donor(donor_tables: dict):
#     """Concatenate tissue-specific candidate peaks such that each donor has a 
#     dataframe with all peaks."""
#     donor_beds = []

#     for donor in donor_tables:
#         beds = []
#         for table in donor_tables[donor]:
#             bed = pd.read_csv(table, sep="\t", 
#                               usecols = ["chrm","start","end","peak","RPM","strand","classification"])
#             bed = bed[bed["classification"].isin(["FP","UNK","SOM_clonal","SOM_private"])]
#             bed = bed[["chrm","start","end","peak","RPM","strand"]]
#             beds.append(bed)
#         donor_bed = pd.concat(beds)
#         donor_bed = donor_bed.sort_values(by=["chrm","start"], ascending=True)
#         donor_beds.append(donor_bed)

#     return donor_beds

def merge_and_filter_peaks(donor_beds: list):
    """Merge peaks across tissues for each donor."""
    merge_beds = []

    for donor_bed in donor_beds:
        donor_bt = BedTool.from_dataframe(donor_bed)
        merge_bt = donor_bt.merge(s=True, c=[4,5,6], o=["distinct","mean","distinct"])
        merge_bed = merge_bt.to_dataframe()
        merge_beds.append(merge_bed)

    return merge_beds   

def intersect_donor_peaks(merge_beds: list):
    """Intersect all peaks from each donor. Noise peaks that are recurrent in
    multiple donors will be output as error-prone regions."""
    x = BedTool()
    merge_fns = [BedTool.from_dataframe(merge_bed).fn for merge_bed in merge_beds]
    regions_bt = x.multi_intersect(i=merge_fns)
    regions_bt_merge = regions_bt.merge(c=4, o="max")
    err_regions = regions_bt_merge.to_dataframe(disable_auto_names=True, header=None)
    # err_regions = err_regions[err_regions[3] > 1] # peak in multiple files
    err_regions = err_regions[[0,1,2,3]]
    err_regions[4] = err_regions[3] / len(merge_beds) # % of donors with peak
    return err_regions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to generate error-prone regions.")
    parser.add_argument("donor_list", help="List of donors to generate regions from")
    parser.add_argument("resdir", help="Parent directory with pipeline results")
    parser.add_argument("outfile", help="Path to output file of error-prone regions")
    args = parser.parse_args()

    with open(args.donor_list, "r") as f:
        donors = [l.strip() for l in f.readlines()]

    donor_beds = []
    for donor in donors:
        peak_file = glob.glob(
            f"{os.path.normpath(args.resdir)}/all*{donor}*classified_peaks.bed")
        
        donor_bed = pd.read_csv(peak_file[0],
                                sep="\t",
                                usecols=[0,1,2,3,4,5,6,7]).sort_values(by=["chrm","start"])
        donor_bed = donor_bed[(donor_bed["classification"] == "UNK") |
                              (donor_bed["classification"].str.contains("SOM")) |
                              (donor_bed["classification"] == "FP")]
        donor_bed = donor_bed[donor_bed["filter_reason"] != "nearby_peak"]
        donor_beds.append(donor_bed)
    
    # donor_beds = concat_bed_by_donor(donor_tables)
    merge_beds = merge_and_filter_peaks(donor_beds)
    err_regions = intersect_donor_peaks(merge_beds)
    err_regions.to_csv(f"{args.outfile}", sep="\t", header=False, index=False)
