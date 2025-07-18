#!/bin/bash

# Calculate precision, recall, and F1 score for each sample.

pdir=$1 # Parent directory containing all pipeline output
truth_name=$2 # Name of ground truth set
truth_peaks=$3 # Ground truth set (GTS) of peaks

for call_peaks in $(ls ${pdir}/*/analysis/*big_table_filter_reasons.tsv); do
    sample=$(basename "${call_peaks}" _big_table_filter_reasons.tsv)
    outdir=$(dirname "${call_peaks}")
    echo $sample

    num_called=$(grep -E 'UNK|SOM' "${call_peaks}" | wc -l)
    num_truth=$(wc -l "${truth_peaks}" | awk '{print $1}')
    both_peaks="${outdir}/${sample}_big_table_filter_${truth_name}_noerror.tsv"

    # Number of true peaks found - don't recount if 1 called peak overlaps multiple true peaks
    tp=$(grep -E 'UNK|SOM' "${both_peaks}" | awk '{print $4}' | sort | uniq | wc -l) 

    fp=$((${num_called} - ${tp}))
    fn=$((${num_truth} - ${tp}))

    precision=$(echo "scale=4; ${tp} / (${tp} + ${fp})" | bc) 
    recall=$(echo "scale=4; ${tp} / (${tp} + ${fn})" | bc) 
    f1=$(echo "scale=4; 2 * (${precision} * ${recall}) / (${precision} + ${recall})" | bc)

    outfile="${outdir}/${sample}_f1_score_${truth_name}_noerror.txt"
    echo -e "${sample}\n" > "${outfile}"
    echo "Number of UNK/SOM peaks called by HAT-seq: ${num_called}" >> "${outfile}"
    echo "Number of peaks in ground truth set: ${num_truth}" >> "${outfile}"
    echo "Number of true peaks found by HAT-seq: ${tp}" >> "${outfile}"
    echo "Precision: ${precision}" >> "${outfile}"
    echo "Recall: ${recall}" >> "${outfile}"
    echo "F1-score: ${f1}" >> "${outfile}"
done
