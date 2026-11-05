import pysam

ibam=pysam.AlignmentFile(snakemake.input.trim_uniq_sorted_bam,'rb')
ibed=open(snakemake.input.peaks, 'r')
obed=open(snakemake.output.unique_start_positions, 'w')

header = ["#peak_name", "num_unique_reads", "num_usp;depth", "num_uep;depth", "max_usp_distance"]
obed.write("\t".join(header) + "\n")

# ibam_total_reads = sum([i[3] for i in ibam.get_index_statistics()])

for peak_line in ibed:
	if peak_line.startswith("#"):
		continue

	peak = peak_line.strip().split("\t")
	peak_iter = ibam.fetch(peak[0], int(peak[1]), int(peak[2]))

	# Unique start positions used as proxy for number of unique PCR templates
	# Unique end positions mark where read ends
	this_peak_usp_dictionary = {}
	this_peak_uep_dictionary = {}

	# this_peak_unique_reads = 0

	for read in peak_iter:
		if peak[5] == "-":
			if read.is_reverse:
				if read.cigartuples[-1][0] == 4:
					skey = read.reference_name + ":" + str(read.reference_end + read.cigartuples[-1][1])
				else: 
					skey = read.reference_name + ":" + str(read.reference_end)

				if not skey in this_peak_usp_dictionary.keys():
					this_peak_usp_dictionary[skey] = 0

				ekey = read.reference_name + ":" + str(read.reference_start)
				if not ekey in this_peak_uep_dictionary.keys():
					this_peak_uep_dictionary[ekey] = 0

				this_peak_usp_dictionary[skey] += 1
				this_peak_uep_dictionary[ekey] += 1
				# this_peak_unique_reads += 1
	
		elif peak[5] == "+":
			if not read.is_reverse:
				if read.cigartuples[0][0] == 4:
					key = read.reference_name + ":" + str(read.reference_start - read.cigartuples[0][1])
				else:
					key = read.reference_name + ":" + str(read.reference_start)
				
				if not key in this_peak_usp_dictionary.keys():
					this_peak_usp_dictionary[key] = 0

				ekey = read.reference_name + ":" + str(read.reference_end)
				if not ekey in this_peak_uep_dictionary.keys():
					this_peak_uep_dictionary[ekey] = 0

				this_peak_usp_dictionary[key] += 1
				this_peak_uep_dictionary[ekey] += 1
				# this_peak_unique_reads += 1
	
	this_peak_unique_reads = sum(
		[this_peak_usp_dictionary[key] for key in this_peak_usp_dictionary.keys()])
	
	skey_list = this_peak_usp_dictionary.keys()
	usp_list = [int(s.split(":")[1]) for s in skey_list]
	left_usp = min(usp_list)
	right_usp = max(usp_list)
	usp_dist = right_usp - left_usp

	sduplicates_list = [this_peak_usp_dictionary[key] for key in this_peak_usp_dictionary.keys()]
	sduplicates_list.sort(reverse=True)

	eduplicates_list = [this_peak_uep_dictionary[key] for key in this_peak_uep_dictionary.keys()]
	eduplicates_list.sort(reverse=True)

	ostring = (peak[3] + "\t" + 
			#    str(ibam_total_reads) + "\t" + 
			#    str(this_peak_total_reads) + "\t" + 
			   str(this_peak_unique_reads) + "\t" +
			   str(len(this_peak_usp_dictionary.keys())) + ";" + 
			   ",".join(map(str, sduplicates_list)) + "\t" +
			   str(len(this_peak_uep_dictionary.keys())) + ";" + 
			   ",".join(map(str, eduplicates_list)) + "\t" +
			   str(usp_dist) + "\n")
	obed.write(ostring) 

obed.close()
ibed.close()
