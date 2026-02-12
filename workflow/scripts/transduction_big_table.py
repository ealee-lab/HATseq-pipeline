import pandas as pd

readIDs = pd.read_table(snakemake.input.peak_readID_list, sep="\t", header=None,
						names=["peak_name","read"])
soft3p = pd.read_table(snakemake.input.softclip3p_trimmed_bed, sep="\t", header=None,
					   names=["3p_soft_clip","start_s3p","end_s3p",
							  "read","mapq_s3p","strand_s3p"])
readIDs = readIDs.merge(soft3p[["read","3p_soft_clip"]], how="left", on="read")
del soft3p 

soft5p = pd.read_table(snakemake.input.softclip5p_trimmed_bed, sep="\t", header=None,
					   names=["5p_soft_clip","start_s5p","end_s5p",
							  "read","mapq_s5p","strand_s5p"])
readIDs = readIDs.merge(soft5p[["read","5p_soft_clip"]], how="left", on="read")
del soft5p 

hard5p = pd.read_table(snakemake.input.hardclip5p_bed, sep="\t", header=None, 
					   names=["chr","start","end","read","null","strand",
							  "chr_ref","start_ref","end_ref","5p_hard_clip","peaks"])
readIDs = readIDs.merge(hard5p[["read","5p_hard_clip"]], how="left", on="read")
del hard5p

hard3p = pd.read_table(snakemake.input.hardclip3p_bed, sep="\t", header=None, 
					   names=["chr","start","end","read","null","strand","chr_ref","start_ref","end_ref","3p_hard_clip","peaks"])
readIDs = readIDs.merge(hard3p[["read","3p_hard_clip"]], how="left", on="read")
del hard3p

readIDs = readIDs.fillna('')
big_table = readIDs.groupby('peak_name').agg({
	'read': lambda x: len(list(set(x))),
	'3p_soft_clip': lambda x: ','.join(list(set(x))),
	'5p_soft_clip': lambda x: ','.join(list(set(x))),
	'5p_hard_clip': lambda x: ','.join(list(set(x))),
	'3p_hard_clip': lambda x: ','.join(list(set(x)))}).reset_index()
del readIDs

custom_ref_peak = pd.read_table(snakemake.input.transduction_bed, sep="\t", 
								usecols=[3,4], names=["custom_ref_peak","peak_name"])
custom_ref_peak["peak_name"] = custom_ref_peak["peak_name"].str.split(",")
custom_ref_peak = custom_ref_peak.explode("peak_name")
big_table = big_table.merge(custom_ref_peak, how="left", on="peak_name")
del custom_ref_peak

big_table = big_table.replace(r',{2,}', '', regex=True)
big_table = big_table.replace(r'^,', '', regex=True)
big_table.to_csv(snakemake.output.transduction_big_table, sep="\t", 
				 index=False, header=True, na_rep="")
