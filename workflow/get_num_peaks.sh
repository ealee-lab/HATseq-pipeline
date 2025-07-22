#!/bin/bash

# Calculate the number of peaks filtered out at each step.

pdir=$1 # Parent directory containing all pipeline output
truth_name=$2 # Name of ground truth set

for both_peaks in $(ls ${pdir}/*/analysis/*big_table_filter_${truth_name}.tsv); do
    sample=$(basename "${both_peaks}" _big_table_filter_${truth_name}.tsv)
    outdir=$(dirname "${both_peaks}")

    echo -e "\nBEFORE ERROR-PRONE FILTER"
    echo "${sample} ${truth_name}"
    awk '{print $4,$7,$8}' "${both_peaks}" | sort | uniq | awk '{print $2,$3}' | sort | uniq -c

    noerror=$(echo "${both_peaks}" | sed 's/.tsv/_noerror.tsv/g')
    echo -e "\nAFTER ERROR-PRONE FILTER"
    awk '{print $4,$7,$8}' "${noerror}" | sort | uniq | awk '{print $2}' | sort | uniq -c | grep -E 'FP|UNK|SOM'
done
