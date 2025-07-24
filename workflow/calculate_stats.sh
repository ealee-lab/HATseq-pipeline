#!/bin/bash

# Calculate precision, recall, and F1 score for each sample.

peak_file=$1 # Path to input peaks intersected with ground truth set (GTS)
noerr_file=$2 # Path to input peaks intersected with GTS and error-prone filtered
out_file=$3 # Path to output file with statistics
true_peaks=$4 # Path to GTS
names=$5 # List of peaks to remove from total count

num_true=$(wc -l "${true_peaks}" | awk '{print $1}')

#### BEFORE ERROR-PRONE FILTER ####

# If names of insertions are given, remove them from the count
# e.g. If using somatic GTS, remove germline GTS insertions from count to subset
if [[ $# == 5 ]]; then

    # Number of called peaks (don't recount if 1 called peak overlaps >1 true peaks)
    num_called=$(grep -v -f "${names}" "${peak_file}" | grep -E 'UNK|SOM' | \
        awk -v FS='\t' '{print $4}' | sort | uniq | wc -l)

    # Number of true peaks found (i.e. number of peaks that overlap with GTS)
    tp=$(grep -v -f "${names}" "${peak_file}" | grep -E 'UNK|SOM' | \
        awk -v FS='\t' -v OFS='\t' '{if ($15 > 0) print $4,$7}' | \
        awk -v FS='\t' '{print $1}' | sort | uniq | wc -l) 

else
    num_called=$(grep -E 'UNK|SOM' "${peak_file}" | \
        awk -v FS='\t' '{print $4}' | sort | uniq | wc -l)

    tp=$(grep -E 'UNK|SOM' "${peak_file}" | \
        awk -v FS='\t' -v OFS='\t' '{if ($15 > 0) print $4,$7}' | \
        awk -v FS='\t' '{print $1}' | sort | uniq | wc -l) 
fi

fp=$((${num_called} - ${tp}))
fn=$((${num_true} - ${tp}))

precision=$(echo "scale=4; ${tp} / (${tp} + ${fp})" | bc) 
recall=$(echo "scale=4; ${tp} / (${tp} + ${fn})" | bc) 
f1=$(echo "scale=4; 2 * (${precision} * ${recall}) / (${precision} + ${recall})" | bc)

echo -e "Number of true insertions: ${num_true}\n" > "${out_file}"
echo -e "BEFORE ERROR-PRONE FILTER: ${sample}" >> "${out_file}"
echo "Number of UNK/SOM peaks called by HAT-seq: ${num_called}" >> "${out_file}"
echo "Number of true peaks found by HAT-seq: ${tp}" >> "${out_file}"
echo "Precision: ${precision}" >> "${out_file}"
echo "Recall: ${recall}" >> "${out_file}"
echo "F1-score: ${f1}" >> "${out_file}"


#### AFTER ERROR-PRONE FILTER ####

if [[ $# == 5 ]]; then
    num_called_noerr=$(grep -v -f "${names}" "${noerr_file}" | grep -E 'UNK|SOM' | \
        awk -v FS='\t' '{print $4}' | sort | uniq | wc -l)

    tp_noerr=$(grep -v -f "${names}" "${noerr_file}" | grep -E 'UNK|SOM' | \
        awk -v FS='\t' -v OFS='\t' '{if ($15 > 0) print $4,$7}' | \
        awk -v FS='\t' '{print $1}' | sort | uniq | wc -l) 
else
    num_called_noerr=$(grep -E 'UNK|SOM' "${noerr_file}" | \
        awk -v FS='\t' '{print $4}' | sort | uniq | wc -l)

    tp_noerr=$(grep -E 'UNK|SOM' "${noerr_file}" | \
        awk -v FS='\t' -v OFS='\t' '{if ($15 > 0) print $4,$7}' | \
        awk -v FS='\t' '{print $1}' | sort | uniq | wc -l) 
fi

fp_noerr=$((${num_called_noerr} - ${tp_noerr}))
fn_noerr=$((${num_true} - ${tp_noerr}))

precision_noerr=$(echo "scale=4; ${tp_noerr} / (${tp_noerr} + ${fp_noerr})" | bc) 
recall_noerr=$(echo "scale=4; ${tp_noerr} / (${tp_noerr} + ${fn_noerr})" | bc) 
f1_noerr=$(echo "scale=4; 2 * (${precision_noerr} * ${recall_noerr}) / (${precision_noerr} + ${recall_noerr})" | bc)

echo -e "\nAFTER ERROR-PRONE FILTER: ${sample}" >> "${out_file}"
echo "Number of UNK/SOM peaks called by HAT-seq: ${num_called_noerr}" >> "${out_file}"
echo "Number of true peaks found by HAT-seq: ${tp_noerr}" >> "${out_file}"
echo "Precision: ${precision_noerr}" >> "${out_file}"
echo "Recall: ${recall_noerr}" >> "${out_file}"
echo "F1-score: ${f1_noerr}" >> "${out_file}"
