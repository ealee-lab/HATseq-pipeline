import pysam

ibam=pysam.AlignmentFile(snakemake.input.peak_sorted_bam,'rb')
ibed=open(snakemake.input.peaks, 'r')

peak_readID_list = []

for peak_line in ibed:
	if peak_line.startswith("#"):
		continue

	peak = peak_line.strip().split("\t")
	peak_iter = ibam.fetch(peak[0], int(peak[1]), int(peak[2]))

	for read in peak_iter:
		if (peak[5] == "-" and read.is_reverse) or (peak[5] == "+" and not read.is_reverse):
			peak_readID_list.append([peak[3], read.query_name])

with open(snakemake.output.peak_readID_list, 'w') as list:
	for peak_read in peak_readID_list:
		list.write("\t".join(peak_read) + "\n")
