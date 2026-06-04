import pysam

def reverse_complement(seq):
	complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
	return ''.join(complement[base] for base in reversed(seq))

otable3p=open(snakemake.output.endpoint3p_table, 'w')

ibam=pysam.AlignmentFile(snakemake.input.peak_sorted_bam, 'rb')
for read in ibam.fetch(until_eof=True):
	if read.is_reverse:
		ostring = (read.query_name + "\t" + 
			 	   read.reference_name + "\t" + 
			 	   str(read.reference_start) + "\t" + 
				   "-")
		
		if "S" in read.cigarstring:
			if read.cigartuples[0][0] == 4: # 3' soft clip
				clip_seq = reverse_complement(
					read.query_sequence[0:read.query_alignment_start])
			elif read.cigartuples[len(read.cigartuples)-1][0] == 4: # 5' soft clip
				clip_seq = ""
				
			ostring = ostring + "\t" + clip_seq

		otable3p.write(ostring + "\n")

	else:
		ostring = (read.query_name + "\t" + 
			 	   read.reference_name + "\t" + 
			 	   str(read.reference_end) + "\t" + 
				   "+")
		
		if "S" in read.cigarstring:
			if read.cigartuples[0][0] == 4: # 5' soft clip
				clip_seq = ""
			elif read.cigartuples[len(read.cigartuples)-1][0] == 4: # 3' soft clip
				clip_seq = read.query_sequence[read.query_alignment_end:read.query_length]

			ostring = ostring + "\t" + clip_seq
		
		otable3p.write(ostring + "\n")

otable3p.close()
