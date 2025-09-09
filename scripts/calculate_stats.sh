#!/bin/bash

# Calculate precision, recall, and F1 score for each sample.

peak_file=$1 # Path to filtered_peaks.tsv intersected with ground truth set (GTS)
true_peaks=$2 # Path to GTS
out_file=$3 # Path to output file with statistics

# Number of false positives (i.e. number of peaks that don't overlap GTS) 
fp=$(grep -E 'UNK|SOM' "${peak_file}" | \
    awk -v FS='\t' '{if ($NF == "-NA-") print $4}' | wc -l)

# Number of true peaks found (i.e. number of peaks that overlap with GTS)
# (Don't recount if 1 called peak overlaps >1 true peaks)
tp=$(awk -v FS='\t' '{print $4}' "${true_peaks}" | grep -f - "${peak_file}" | \
    grep -E 'UNK|SOM' | awk -v FS='\t' '{if ($NF != "-NA-") print $4}' | \
    sort | uniq | wc -l) 

num_called=$((${tp} + ${fp}))
num_true=$(wc -l "${true_peaks}" | awk '{print $1}')
fn=$((${num_true} - ${tp}))

precision=$(echo "scale=4; ${tp} / (${tp} + ${fp})" | bc) 
recall=$(echo "scale=4; ${tp} / (${tp} + ${fn})" | bc) 
f1=$(echo "scale=4; 2 * (${precision} * ${recall}) / (${precision} + ${recall})" | bc)

echo -e "Number of true insertions: ${num_true}" > "${out_file}"
echo "Number of UNK/SOM peaks called by HAT-seq: ${num_called}" >> "${out_file}"
echo "Number of true peaks found by HAT-seq: ${tp}" >> "${out_file}"
echo "Precision: ${precision}" >> "${out_file}"
echo "Recall: ${recall}" >> "${out_file}"
echo "F1-score: ${f1}" >> "${out_file}"
