import pandas as pd
import statistics


def calculate_breakpoint_concordance(x):
	bps = [int(y) for y in x if y >= 0]

	if len(bps) > 0:
		bp_near = [d for d in bps if d <= 8]
		bp_pct = len(bp_near) / len(bps)
	else:
		bp_pct = -1
	return bp_pct


peaks = pd.read_table(
	snakemake.input.peaks, sep="\t", skiprows=1, header=None, 
	names=["chr","start","end","peak_name","shape","strand"])

readIDs = pd.read_table(
	snakemake.input.peak_readID_list, sep="\t", header=None, names=["peak_name","read"])
peaks = peaks.merge(readIDs, how="left", on="peak_name")
del readIDs

# Ntag = pd.read_table(snakemake.input.Ntag_list, sep="\t", 
# 					 header=None, names=["read","Ntag"])
# peaks = peaks.merge(Ntag, how="left", on="read").fillna("")
# del Ntag

gmotif = pd.read_table(snakemake.input.pass_gmotif_list, sep="\t", header=None, names=["read"])
gmotif["gmotif"] = "pass_" + gmotif["read"]
peaks = peaks.merge(gmotif, how="left", on="read").fillna("fail")
del gmotif

endpoints = pd.read_table(
	snakemake.input.endpoint3p_table, sep="\t", header=None, 
	names=["read","chr","endpoint","strand","clip_seq"])
peaks = peaks.merge(endpoints, how="left", on=["read","chr","strand"]).fillna(".")
del endpoints

# Delete rows where clipped read is incorrectly matched to peak
peaks = peaks[(peaks["endpoint"] >= peaks["start"]) & (peaks["endpoint"] <= peaks["end"])]

polyA = pd.read_table(
	snakemake.input.polyT_reads, sep="\t", header=None, 
	usecols=[0,1,2,3], names=["read","chr","endpoint","strand"])
polyA["polyA"] = "pass_" + polyA['read']
peaks = peaks.merge(polyA, how="left", on=["read","chr","endpoint","strand"]).fillna("fail")
del polyA

breakpoints = peaks[peaks["polyA"] != "fail"].groupby("peak_name").agg(
	breakpoint=('endpoint', lambda x: statistics.mode(x))
).reset_index()
peaks = peaks.merge(breakpoints, how="left", on="peak_name").fillna(-1)

# Calculate distance from breakpoint for select set of reads
peaks["bp_dist"] = -1

# Distance for all clipped reads
peaks.loc[
	(peaks["breakpoint"] > 0) &
	(peaks["clip_seq"] != "."), "bp_dist"] = abs(peaks["endpoint"] - peaks["breakpoint"])

# Distance for all reads extending past the breakpoint
peaks.loc[
	(peaks["breakpoint"] > 0) & 
	(peaks["endpoint"] >= peaks["breakpoint"]) &
	(peaks["strand"] == "+"), "bp_dist"] = abs(peaks["endpoint"] - peaks["breakpoint"])
peaks.loc[
	(peaks["breakpoint"] > 0) & 
	(peaks["endpoint"] <= peaks["breakpoint"]) &
	(peaks["strand"] == "-"), "bp_dist"] = abs(peaks["endpoint"] - peaks["breakpoint"])

big_table = peaks.groupby(["chr","start","end","peak_name","shape","strand"]).agg({
	'read': lambda x: len(list(set(x))),
	# 'Ntag': lambda x: ",".join([y for y in list(set(x)) if (len(y) == 2 or len(y) == 4 or len(y) == 6)]),
	'gmotif': lambda x: len(list(set([y for y in x if y != "fail"]))),
	'polyA': lambda x: len(list(set([y for y in x if y != "fail"]))),
	'bp_dist': lambda x: calculate_breakpoint_concordance(x)
}).reset_index()
del peaks

max_depth_distance = pd.read_table(snakemake.input.max_depth,sep="\t", 
								   header=None, names=['peak_name','distance'])
big_table = big_table.merge(max_depth_distance, how="left", on="peak_name")
del max_depth_distance

nearby_peak = pd.read_table(snakemake.input.nearby_peaks, sep="\t", 
							names=['nearest_peak','peak_name'])
nearby_peak['peak_name'] = nearby_peak['peak_name'].str.split(pat="=")
nearby_peak = nearby_peak.explode('peak_name')
big_table = big_table.merge(nearby_peak, how="left", on="peak_name").fillna(".")
del nearby_peak

satellites = pd.read_table(snakemake.input.satellite_intersect, sep="\t", 
						   names=['peak_name','Satellite'])
big_table = big_table.merge(satellites, how="left", on="peak_name").fillna(".")
del satellites

segdups = pd.read_table(snakemake.input.segdup_intersect, sep="\t", 
						names=['peak_name','SegDup'])
big_table = big_table.merge(segdups, how="left", on="peak_name").fillna(".")
del segdups

homopolymers = pd.read_table(snakemake.input.homopolymer_intersect, sep="\t", 
							 names=['peak_name','homopolymer'])
big_table = big_table.merge(homopolymers, how="left", on="peak_name").fillna(".")
del homopolymers

intersect = pd.read_table(snakemake.input.intersect_annotated, sep="\t", header=0)
big_table = big_table.merge(
	intersect[['peak_name','RepeatMasker','Evrony_KR',
			   '1000_Genomes_Project','gnomAD','NyuWa','xTea',
			   'HGSVC3','HGSVC3-MELT-LRA','1019_ONT']], 
	how="left", on="peak_name")
del intersect

usp = pd.read_table(snakemake.input.unique_start_positions, sep="\t", header=0)
usp.columns = usp.columns.str.strip("#")
usp["RPM"] = (usp['num_peak_reads'] / (usp['num_bam_reads'])) * 1000000
usp["unique_read_ratio"] = usp["num_peak_unique_reads"] / usp["num_peak_reads"]
big_table = big_table.merge(usp, how="left", on="peak_name")
del usp

big_table = big_table.replace(r',{2,}', '', regex=True)
big_table.to_csv(snakemake.output.big_table, sep="\t", index=False, header=True, na_rep="")
