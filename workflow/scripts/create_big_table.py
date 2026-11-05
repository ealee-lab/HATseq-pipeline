import sys
import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w", buffering=1)

peaks = pd.read_table(
	snakemake.input.peaks, sep="\t", skiprows=1, header=None, 
	names=["chrm","start","end","peak_name","shape","strand"])
peaks = peaks.astype(
	{"chrm": str, "start": int, "end": int, "peak_name": str, "shape": str, "strand": str})

readIDs = pd.read_table(
	snakemake.input.peak_readID_list, sep="\t", header=None, names=["peak_name","read"])
peaks = peaks.merge(readIDs, how="left", on="peak_name")
del readIDs

num_bam_reads = len(peaks)

gmotif = pd.read_table(snakemake.input.pass_gmotif_list, sep="\t", header=None, names=["read"])
gmotif["gmotif"] = "pass_" + gmotif["read"]
peaks = peaks.merge(gmotif, how="left", on="read").fillna("fail")
del gmotif

polyA = pd.read_table(snakemake.input.polyT_reads, sep="\t", header=None,
					  usecols=[0,1,2,3], names=["read","chrm","clippoint","strand"])
uniq = pd.read_table(snakemake.input.uniq_readID_list, sep="\t", header=None, names=["read"])
polyA = polyA.merge(uniq, how="inner", on="read")
polyA["polyA"] = "pass_" + polyA['read']
peaks = peaks.merge(polyA, how="left", on=["read","chrm","strand"])
del polyA

peaks["polyA"].fillna("fail", inplace=True)
peaks.loc[peaks["clippoint"].isna(), "clippoint"] = peaks["end"] # arbitrary for reads without clip
peaks["clippoint"] = peaks["clippoint"].astype(int)

# Delete rows where polyA read is incorrectly matched to peak
peaks = peaks[(peaks["clippoint"] >= peaks["start"]) & (peaks["clippoint"] <= peaks["end"])]

big_table = peaks.groupby(["chrm","start","end","peak_name","shape","strand"]).agg(
	num_reads=('read', lambda x: len(list(set(x)))),
	num_gmotif_reads=('gmotif', lambda x: len(list(set([y for y in x if y != "fail"])))),
	num_polyA_reads=('polyA', lambda x: len(list(set([y for y in x if y != "fail"]))))
	# 'bp_dist': lambda x: calculate_breakpoint_concordance(x)
).reset_index()
del peaks

big_table["peak_width"] = big_table["end"] - big_table["start"]
big_table["gmotif_percent"] = (big_table["num_gmotif_reads"] / big_table["num_reads"]) * 100
big_table["polyA_percent"] = (big_table["num_polyA_reads"] / big_table["num_reads"]) * 100

usp = pd.read_table(snakemake.input.unique_start_positions, sep="\t", header=0)
usp.columns = usp.columns.str.strip("#")
big_table = big_table.merge(usp, how="left", on="peak_name")
del usp

big_table["num_templates"] = big_table["num_usp;depth"].str.split(";").str[0].astype(int)
big_table["num_endpoints"] = big_table["num_uep;depth"].str.split(";").str[0].astype(int)

template1 = big_table["num_usp;depth"].str.split(";").str[1].str.split(",").str[0]
template2 = big_table["num_usp;depth"].str.split(";").str[1].str.split(",").str[1]
template2.loc[template2.isna()] = template1 # template_ratio = 1 if only 1 template
big_table["template_ratio"] = template2.astype(int) / template1.astype(int)

big_table["start_end_ratio"] = big_table["num_templates"] / big_table["num_endpoints"]
big_table["unique_read_ratio"] = big_table["num_unique_reads"] / big_table["num_reads"]
big_table["RPM"] = (big_table["num_reads"] / num_bam_reads) * 1000000
big_table["TPM"] = (big_table["num_templates"] / num_bam_reads) * 1000000

big_table.drop(columns=["num_usp;depth","num_uep;depth"], inplace=True)

# max_depth_distance = pd.read_table(snakemake.input.max_depth,sep="\t", 
# 								   header=None, names=['peak_name','distance'])
# big_table = big_table.merge(max_depth_distance, how="left", on="peak_name")
# del max_depth_distance

nearby_peak = pd.read_table(
	snakemake.input.nearby_peaks, sep="\t", names=['nearest_peak','peak_name'])
nearby_peak['peak_name'] = nearby_peak['peak_name'].str.split(pat="=")
nearby_peak = nearby_peak.explode('peak_name')
big_table = big_table.merge(nearby_peak, how="left", on="peak_name").fillna(".")
del nearby_peak

big_table = big_table.merge(big_table[['peak_name','RPM']], how="left", 
							left_on="nearest_peak", right_on="peak_name", suffixes=("", "_nearest"))
big_table["nearest_RPM_ratio"] = big_table["RPM"] / big_table["RPM_nearest"]
big_table["nearest_RPM_ratio"].fillna(1, inplace=True)
big_table.drop(columns="peak_name_nearest", inplace=True)

satellites = pd.read_table(snakemake.input.satellite_intersect, sep="\t", 
						   names=['peak_name','satellite'])
big_table = big_table.merge(satellites, how="left", on="peak_name").fillna(".")
del satellites

segdups = pd.read_table(snakemake.input.segdup_intersect, sep="\t", 
						names=['peak_name','segdup'])
big_table = big_table.merge(segdups, how="left", on="peak_name").fillna(".")
del segdups

homopolymers = pd.read_table(snakemake.input.homopolymer_intersect, sep="\t", 
							 names=['peak_name','homopolymer'])
big_table = big_table.merge(homopolymers, how="left", on="peak_name").fillna(".")
del homopolymers

intersect = pd.read_table(snakemake.input.intersect_annotated, sep="\t", header=0).fillna(".")
big_table = big_table.merge(
	intersect[['peak_name','nearest_KR','nearest_KNR']], 
	how="left", on="peak_name")
del intersect

big_table["nearest_KR_dist"] = big_table["nearest_KR"].str.split(",").str[0].str.split(";").str[1]
big_table.loc[big_table["nearest_KR"] == ".", "nearest_KR_dist"] = 999
big_table["nearest_KNR_dist"] = big_table["nearest_KNR"].str.split(",").str[0].str.split(";").str[1]
big_table.loc[big_table["nearest_KNR"] == ".", "nearest_KNR_dist"] = 999

big_table = big_table.replace(r',{2,}', '', regex=True)
big_table.to_csv(snakemake.output.big_table, sep="\t", index=False, header=True, na_rep=".")
