import os
import glob
import argparse
import pandas as pd
import pybedtools
from pybedtools import BedTool

# def concat_and_merge_annotations(kr: list, knr: list):
#     """Concatenate and merge all KR and KNR annotations, generating one dataframe
#     for each. Annotation files should be in BED format, with KR annotations 
#     specifically in BED6 format."""
#     all_kr = pd.read_csv(kr[0], sep="\t", usecols=[0,1,2,3,4,5],
#                          names=["chr","start","end","name","score","strand"])
#     for ref in kr:
#         kr_df = pd.read_csv(ref, sep="\t", usecols=[0,1,2,3,4,5], header=None)
#         all_kr = pd.concat([all_kr, kr_df])

#     all_knr = pd.read_csv(knr[0], sep="\t", usecols=[0,1,2], 
#                           names=["chr","start","end"])
#     for ann in knr:
#         knr_df = pd.read_csv(ann, sep="\t", usecols=[0,1,2], header=None)
#         all_knr = pd.concat([all_knr, knr_df])

#     all_kr = all_kr.sort_values(by=["chr","start"])
#     all_knr = all_knr.sort_values(by=["chr","start"])

#     kr_bt = BedTool.from_dataframe(all_kr)
#     knr_bt = BedTool.from_dataframe(all_knr)

#     kr_merge_bt = kr_bt.merge(s=True, c=[4,5,6], o=["distinct","mean","distinct"])
#     knr_merge_bt = knr_bt.merge()

#     kr_merge = kr_merge_bt.to_dataframe()
#     knr_merge = knr_merge_bt.to_dataframe()
#     return kr_merge, knr_merge

def concat_bed_by_donor(donor_tables: dict):
    """Concatenate tissue-specific UNK/SOM peaks such that each donor has a 
    dataframe with all peaks."""
    donor_beds = []

    for donor in donor_tables:
        beds = []
        for table in donor_tables[donor]:
            bed = pd.read_csv(table, sep="\t", 
                              usecols = ["chrm","start","end","peak","RPM","strand","classification"])
            bed = bed[(bed["classification"] == "UNK") | (bed["classification"].str.contains("SOM"))]
            bed = bed[["chrm","start","end","peak","RPM","strand"]]
            beds.append(bed)
        donor_bed = pd.concat(beds)
        donor_bed = donor_bed.sort_values(by=["chrm","start"], ascending=True)
        donor_beds.append(donor_bed)

    return donor_beds

def merge_and_filter_peaks(donor_beds: list):
    """Merge adjacent peaks for each donor. Also, only keep noise peaks
    (peaks with RPM < 50)."""
    merge_beds = []

    for donor_bed in donor_beds:
        donor_bt = BedTool.from_dataframe(donor_bed)
        merge_bt = donor_bt.merge(s=True, c=[4,5,6], o=["distinct","mean","distinct"])
        merge_bed = merge_bt.to_dataframe()
        merge_bed = merge_bed[merge_bed["score"] < 50] # count RPM < 50 as noise
        merge_beds.append(merge_bed)

    return merge_beds   

def intersect_donor_peaks(merge_beds: list):
    """Intersect all peaks from each donor. Noise peaks that are recurrent in
    multiple donors will be output as error-prone regions."""
    x = BedTool()
    merge_fns = [BedTool.from_dataframe(merge_bed).fn for merge_bed in merge_beds]
    regions_bt = x.multi_intersect(i=merge_fns)
    err_regions = regions_bt.to_dataframe(disable_auto_names=True, header=None)
    err_regions = err_regions[err_regions[3] > 1] # peak in multiple files
    err_regions = err_regions[[0,1,2,3]]
    err_regions[4] = err_regions[3] / len(merge_beds) # % of donors with peak
    return err_regions

# def remove_KR_KNR_peaks(err_regions: pd.DataFrame, kr: pd.DataFrame, knr: pd.DataFrame):
#     regions_bt = BedTool.from_dataframe(err_regions)
#     kr_bt = BedTool.from_dataframe(kr)
#     knr_bt = BedTool.from_dataframe(knr)

#     # Intersect -v
#     return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script to generate error-prone regions.")
    parser.add_argument("donor_list", help="List of donors to generate regions from")
    parser.add_argument("directory", help="Parent directory with pipeline outputs")
    parser.add_argument("reference", help="Reference directory with KR/KNR annotations")
    args = parser.parse_args()

    with open(args.donor_list, "r") as f:
        donors = [l.strip() for l in f.readlines()]
    
    donor_tables = {}
    for donor in donors:
        tables = glob.glob(
            f"{os.path.normpath(args.directory)}/*{donor}*/*{donor}*big_table_filter_reasons.tsv")
        
        if tables != []: # if files exist for donor
            donor_tables[donor] = tables
        else:
            continue 

    donor_beds = concat_bed_by_donor(donor_tables)
    merge_beds = merge_and_filter_peaks(donor_beds)
    err_regions = intersect_donor_peaks(merge_beds)
    err_regions.to_csv(f"{args.reference}/ErrorProne/bulk.NIH_Aging.error-prone.bed", sep="\t", header=False, index=False)

    # kr_files = glob.glob(f"{os.path.normpath(args.reference)}/RepeatMasker/*.bed")

    # knr_anns = ["1kgp", "1019_ONT", "gnomAD-SV", "HGSVC3", "nyuwa", "xTea"] # TODO: don't hard code this
    # knr_files = []
    # for ann in knr_anns:
    #     knr_file = glob.glob(f"{os.path.normpath(args.reference)}/{ann}/*.bed")
    #     knr_files += knr_file
            
    # kr_df, knr_df = concat_and_merge_annotations(kr_files, knr_files)