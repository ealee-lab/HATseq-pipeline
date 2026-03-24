import gzip
import pysam

ofasta5p=gzip.open(snakemake.output.softclip5p_fasta, 'wt')
otable5p=open(snakemake.output.hardclip5p_table, 'w')
ofasta3p=gzip.open(snakemake.output.softclip3p_fasta, 'wt')
otable3p=open(snakemake.output.hardclip3p_table, 'w')

def reverse_complement(seq):
	complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
	return ''.join(complement[base] for base in reversed(seq))

ibam=pysam.AlignmentFile(snakemake.input.peak_sorted_bam, 'rb')
for read in ibam.fetch(until_eof=True):
	reverse = (read.flag & 16 == 16) 

	if "S" in read.cigarstring:
		if read.cigartuples[0][0] == 4:
			if reverse:
				ofasta3p.write(">" + read.query_name + "\n")
				ofasta3p.write(
					reverse_complement(read.query_sequence[0:read.query_alignment_start]) + "\n")
			else:
				ofasta5p.write(">" + read.query_name + "\n")
				ofasta5p.write(
					read.query_sequence[0:read.query_alignment_start] + "\n")
		elif read.cigartuples[len(read.cigartuples)-1][0] == 4:
			if reverse:
				ofasta5p.write(">" + read.query_name + "\n")
				ofasta5p.write(
					reverse_complement(read.query_sequence[read.query_alignment_end:read.query_length]) + "\n")
			else:
				ofasta3p.write(">" + read.query_name + "\n")
				ofasta3p.write(
					read.query_sequence[read.query_alignment_end:read.query_length] + "\n")

	if "H" in read.cigarstring:
		if read.cigartuples[0][0] == 5:
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
						otable3p.write(ostring)
					else:
						otable5p.write(ostring)

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
						otable5p.write(ostring)
					else:
						otable3p.write(ostring)

ofasta5p.close()
otable5p.close()
ofasta3p.close()
otable3p.close()
