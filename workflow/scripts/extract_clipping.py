import gzip
import pysam
import pandas as pd

sfasta5p=gzip.open(snakemake.output.softclip5p_fasta, 'wt')
htable5p=open(snakemake.output.hardclip5p_table, 'w')
sfasta3p=gzip.open(snakemake.output.softclip3p_fasta, 'wt')
htable3p=open(snakemake.output.hardclip3p_table, 'w')
stable3p=open(snakemake.output.softclip3p_table, 'w')

def reverse_complement(seq):
	complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
	return ''.join(complement[base] for base in reversed(seq))

names = []
chrs = []
endpoints = []
strands = []
seqs = []

ibam=pysam.AlignmentFile(snakemake.input.peak_sorted_bam, 'rb')
for read in ibam.fetch(until_eof=True):
	reverse = (read.flag & 16 == 16)

	if "S" in read.cigarstring:
		if read.cigartuples[0][0] == 4: # alignment starts with soft clip
			if reverse:
				clip_seq = reverse_complement(
					read.query_sequence[0:read.query_alignment_start])
				sfasta3p.write(">" + read.query_name + "\n" + clip_seq + "\n")
				
				stable3p.write(read.query_name + "\t" + read.reference_name + "\t" + 
				   str(read.reference_start) + "\t-\t" + clip_seq + "\n")
				
				# names.append(read.query_name)
				# chrs.append(read.reference_name)
				# endpoints.append(read.reference_start)
				# strands.append("-")
				# seqs.append(clip_seq)
			else:
				clip_seq = read.query_sequence[0:read.query_alignment_start]
				sfasta5p.write(">" + read.query_name + "\n" + clip_seq + "\n")
				
		elif read.cigartuples[len(read.cigartuples)-1][0] == 4: # alignment ends with soft clip
			if reverse:
				clip_seq = reverse_complement(
					read.query_sequence[read.query_alignment_end:read.query_length])
				sfasta5p.write(">" + read.query_name + "\n" + clip_seq + "\n")
				
			else:
				clip_seq = read.query_sequence[read.query_alignment_end:read.query_length]
				sfasta3p.write(">" + read.query_name + "\n" + clip_seq + "\n")
				
				stable3p.write(read.query_name + "\t" + read.reference_name + "\t" + 
				   str(read.reference_end) + "\t+\t" + clip_seq + "\n")
				
				# names.append(read.query_name)
				# chrs.append(read.reference_name)
				# endpoints.append(read.reference_end)
				# strands.append("+")
				# seqs.append(clip_seq)

	if "H" in read.cigarstring:
		if read.cigartuples[0][0] == 5:
			num_alt_mappings = len(read.get_tag("SA").strip(";").split(";"))
			
			if num_alt_mappings == 1:
				uniq_supp_map = read.get_tag("SA").strip(";").split(";")[0].split(",")
				# print(uniq_supp_map)
				if int(uniq_supp_map[4]) > 20 and int(uniq_supp_map[5]) < 4: 
					ostring = (read.query_name + "\t" + 
							   str(uniq_supp_map[0]) + ":" +  
							   str(uniq_supp_map[1]) + "-" + 
							   str(int(uniq_supp_map[1]) + int(read.reference_length)) + ":" + 
							   uniq_supp_map[2] + "\n")
					if reverse:
						htable3p.write(ostring)
					else:
						htable5p.write(ostring)

		elif read.cigartuples[len(read.cigartuples)-1][0] == 5:
			num_alt_mappings = len(read.get_tag("SA").strip(";").split(";"))
			
			if num_alt_mappings == 1:
				uniq_supp_map = read.get_tag("SA").strip(";").split(";")[0].split(",")
				
				if int(uniq_supp_map[4]) > 20 and int(uniq_supp_map[5]) < 4:
					ostring = (read.query_name + "\t" + 
							   str(uniq_supp_map[0]) + ":" +  
							   str(uniq_supp_map[1]) + "-" + 
							   str(int(uniq_supp_map[1]) + int(read.reference_length)) + ":" + 
							   uniq_supp_map[2] + "\n")
					if reverse:
						htable5p.write(ostring)
					else:
						htable3p.write(ostring)

sfasta5p.close()
htable5p.close()
sfasta3p.close()
htable3p.close()

# df_3p = pd.DataFrame({'read': names, 
# 					  'chr': chrs, 
# 					  'endpoint': endpoints, 
# 					  'strand': strands,
# 					  "sequence": seqs}).sort_values(by=["chr","endpoint"])
# df_3p.to_csv(snakemake.output.softclip3p_table, sep="\t", index=False, header=False)
