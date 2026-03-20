import pandas as pd

Ntag = pd.read_table(snakemake.input.Ntag_list, sep="\t", 
					 header=None, names=["read","Ntag"])
readIDs = pd.read_table(snakemake.input.peak_readID_list, sep="\t", 
						header=None, names=["peak_name","read"])
readIDs = readIDs.merge(Ntag, how="left", on="read")
del Ntag

polyA = pd.read_table(snakemake.input.polyT_reads, sep="\t",header=None, names=["read"])
polyA["polyA"] = "pass_" + polyA['read']
readIDs = readIDs.merge(polyA, how="left", on="read").fillna('fail')
del polyA

gmotif = pd.read_table(snakemake.input.pass_gmotif_list, sep="\t", header=None, names=["read"])
gmotif["gmotif"] = "pass_" + gmotif['read']
readIDs = readIDs.merge(gmotif, how="left", on="read").fillna('fail')
del gmotif

readIDs = readIDs.fillna('')
big_table = readIDs.groupby('peak_name').agg({
	'read': lambda x: len(list(set(x))),
	'Ntag': lambda x: ([y for y in list(set(x)) if (len(y) == 2 or len(y) == 4 or len(y) == 6)]),
	'polyA': lambda x: len(list(set([y for y in x if y != "fail"]))),
	'gmotif': lambda x: len(list(set([y for y in x if y != "fail"])))}).reset_index()
del readIDs

peaks = pd.read_table(snakemake.input.peaks, sep="\t", header=0)
big_table = peaks.merge(big_table, how="left", on="peak_name")
del peaks

chimera = pd.read_table(snakemake.input.chimera, sep="\t", 
						header=None, names=["peak_name","bp_chimera_len","bp_chimera_ratio"])
big_table = big_table.merge(chimera, how="left", on="peak_name")
del chimera

intersect = pd.read_table(snakemake.input.intersect_annotated, sep="\t", header=0)
big_table = big_table.merge(
	intersect[['peak_name','RepeatMasker','Evrony_KR','Homopolymers',
			   '1000_Genomes_Project','gnomAD','NyuWa','xTea',
			   'HGSVC3','HGSVC3-MELT-LRA','1019_ONT']], 
	how="left", on="peak_name")
del intersect  

usp = pd.read_table(snakemake.input.unique_start_positions, sep="\t", header=0)
usp.columns = usp.columns.str.strip("#")
big_table = big_table.merge(usp, how="left", on="peak_name")
del usp

big_table["RPM"] = (big_table['num_peak_reads'] / (big_table['num_bam_reads'])) * 1000000
nearby_peak = pd.read_table(snakemake.input.nearby_peaks, sep="\t", 
							names=['nearest_peak','peak_name'])
nearby_peak['peak_name'] = nearby_peak['peak_name'].str.split(pat="=")
nearby_peak = nearby_peak.explode('peak_name')
big_table = big_table.merge(nearby_peak, how="left", on="peak_name")
del nearby_peak

segdups = pd.read_table(snakemake.input.segdup_intersect, sep="\t", 
						names=['chr','start','end','peak_name','SegDup'])
big_table = big_table.merge(segdups[['peak_name','SegDup']], how="left", on="peak_name")
del segdups

max_depth_distance = pd.read_table(snakemake.input.max_depth,sep="\t", 
								   header=None, names=['peak_name','distance'])
big_table = big_table.merge(max_depth_distance, how="left", on="peak_name")
del max_depth_distance 

big_table = big_table.replace(r',{2,}', '', regex=True)
big_table.to_csv(snakemake.output.big_table, sep="\t", index=False, header=True, na_rep="")
