import pysam
import statistics

ibam=pysam.AlignmentFile(snakemake.input.peak_sorted_bam,'rb')

ibed=open(snakemake.input.peaks, 'r')
obed=open(snakemake.output.unique_start_positions, 'w')
header = ["#peak_name", "num_bam_reads", "num_peak_reads", "num_usp;depth", 
		  "max_neighbor_difference", "median_neighbor_difference", "mean_neighbor_difference"]
obed.write("\t".join(header) + "\n")

ibam_total_reads = sum([i[3] for i in ibam.get_index_statistics()])
peak_readID_list = []

for peak_line in ibed:
	if peak_line.startswith("#"):
		continue

	peak = peak_line.strip().split("\t")
	peak_iter = ibam.fetch(peak[0], int(peak[1]), int(peak[2]))
	this_peak_usp_dictionary = {}

	for read in peak_iter:
		if read.is_reverse:
			key = read.reference_name + ":" + str(read.reference_end)
		else:
			key = read.reference_name + ":" + str(read.reference_start)

		if peak[5] == "-":
			if read.is_reverse:
				if not key in this_peak_usp_dictionary.keys():
					this_peak_usp_dictionary[key] = (0,0)
				peak_readID_list.append([peak[3], read.query_name])

				this_peak_usp_dictionary[key] = (this_peak_usp_dictionary[key][0]+1,
												 this_peak_usp_dictionary[key][1]+1)
		elif peak[5] == "+":
			if not read.is_reverse:
				if not key in this_peak_usp_dictionary.keys():
					this_peak_usp_dictionary[key] = (0,0)
				peak_readID_list.append([peak[3], read.query_name])

				this_peak_usp_dictionary[key] = (this_peak_usp_dictionary[key][0]+1,
												 this_peak_usp_dictionary[key][1]+1)
	
	this_peak_total_reads = sum(
		[this_peak_usp_dictionary[key][1] for key in this_peak_usp_dictionary.keys()])
	duplicates_list = [this_peak_usp_dictionary[key][1] for key in this_peak_usp_dictionary.keys()]
	difference_list = []

	if(len(duplicates_list) > 1):
		for i in range(0,len(duplicates_list)-1):
			difference_list.append(abs(duplicates_list[i] - duplicates_list[i+1]))
	else:
		difference_list = [-1]

	duplicates_list.sort(reverse=True)

	ostring = (peak[3] + "\t" + 
			      str(ibam_total_reads) + "\t" + 
			      str(this_peak_total_reads) + "\t" + 
				  str(len(this_peak_usp_dictionary.keys())) + ";" + 
				  ",".join(map(str, duplicates_list)) + "\t" + 
				  str(max(difference_list)) + "\t" + 
				  str(statistics.median(difference_list)) + "\t" + 
				  str(statistics.mean(difference_list)) + "\n")
	obed.write(ostring) 

with open(snakemake.output.peak_readID_list, 'w') as list:
	for peak_read in peak_readID_list:
		list.write("\t".join(peak_read) + "\n")

obed.close()
ibed.close()
