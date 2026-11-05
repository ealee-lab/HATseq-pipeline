#!/bin/bash

# Calculate precision, recall, and F1 score for each sample.

peak_file=$1 # Path to all_donor_classified_peaks.bed
true_peaks=$2 # Path to GTS
out_file=$3 # Path to output file with statistics
# sample_list=$4 # List of samples

num_true=$(wc -l "${true_peaks}" | awk '{print $1}')
echo -e "Number of true insertions: ${num_true}" > "${out_file}"

samples=$(awk 'NR > 1 \
			{sub(/-(minus|plus)-peak-[0-9]+$/, "", $4); print $4}' "${peak_file}" | \
		  sort | uniq)
# readarray -t samples < "${sample_list}"

for sample in ${samples[@]}; do
	# Number of false positives (i.e. number of peaks that don't overlap GTS) 
	fp=$(grep -w "${sample}" "${peak_file}" | \
		grep -E 'UNK|SOM' | \
		awk -v FS='\t' '{if ($NF == "-NA-") print $4}' | \
		wc -l)

	# Number of true peaks found (i.e. number of peaks that overlap with GTS)
	# Don't re-count if 1 true peak overlaps >1 called peaks
	# Don't re-count if 1 called peak overlaps >1 true peaks
	tp=$(awk -v FS='\t' '{print $4}' "${true_peaks}" | \
		grep -w -f - "${peak_file}" | \
		grep -w "${sample}" | \
		grep -E 'UNK|SOM' | \
		awk '!seen[$NF]++' | \
		awk '!seen[$4]++' | \
		wc -l)

	num_called=$((${tp} + ${fp}))
	fn=$((${num_true} - ${tp}))

	precision=$(echo "scale=4; ${tp} / (${tp} + ${fp})" | bc) 
	recall=$(echo "scale=4; ${tp} / (${tp} + ${fn})" | bc) 
	f1=$(echo "scale=4; 2 * (${precision} * ${recall}) / (${precision} + ${recall})" | bc)

	mkdir -p $(dirname ${out_file})

	echo -e "\n${sample}" >> "${out_file}"
	echo "Number of UNK/SOM peaks called by HAT-seq: ${num_called}" >> "${out_file}"
	echo "Number of true peaks found by HAT-seq: ${tp}" >> "${out_file}"
	echo "Precision: ${precision}" >> "${out_file}"
	echo "Recall: ${recall}" >> "${out_file}"
	echo "F1-score: ${f1}" >> "${out_file}"
done
