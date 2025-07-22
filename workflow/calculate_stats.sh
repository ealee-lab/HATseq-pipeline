#!/bin/bash

# Calculate precision, recall, and F1 score for each sample.

pdir=$1 # Parent directory containing all pipeline output
truth_name=$2 # Name of ground truth set
truth_peaks=$3 # Ground truth set (GTS) of peaks

num_truth=$(wc -l "${truth_peaks}" | awk '{print $1}')

for path in $(ls ${pdir}/); do
    sample=$(basename "${path}")
    outdir=$(realpath "${pdir}")/${sample}/analysis/post
    echo $sample

    # BEFORE ERROR-PRONE FILTER
    call_peaks="${pdir}/${sample}/analysis/${sample}_filtered_peaks.tsv"
    both_peaks="${outdir}/${sample}_big_table_filter_${truth_name}.tsv"

    # Number of called peaks
    num_called=$(grep -E 'UNK|SOM' "${call_peaks}" | wc -l)

    # Number of true peaks found - don't recount if 1 called peak overlaps multiple true peaks
    tp=$(grep -E 'UNK|SOM' "${both_peaks}" | awk '{print $4}' | sort | uniq | wc -l) 

    fp=$((${num_called} - ${tp}))
    fn=$((${num_truth} - ${tp}))

    precision=$(echo "scale=4; ${tp} / (${tp} + ${fp})" | bc) 
    recall=$(echo "scale=4; ${tp} / (${tp} + ${fn})" | bc) 
    f1=$(echo "scale=4; 2 * (${precision} * ${recall}) / (${precision} + ${recall})" | bc)

    outfile="${outdir}/${sample}_f1_score_${truth_name}.txt"
    echo -e "Number of true insertions: ${num_truth}\n" > "${outfile}"
    echo -e "BEFORE ERROR-PRONE FILTER: ${sample}" >> "${outfile}"
    echo "Number of UNK/SOM peaks called by HAT-seq: ${num_called}" >> "${outfile}"
    echo "Number of true peaks found by HAT-seq: ${tp}" >> "${outfile}"
    echo "Precision: ${precision}" >> "${outfile}"
    echo "Recall: ${recall}" >> "${outfile}"
    echo "F1-score: ${f1}" >> "${outfile}"


    # AFTER ERROR-PRONE FILTER
    call_peaks_noerr="${outdir}/${sample}_big_table_filter_noerror.tsv"
    both_peaks_noerr="${outdir}/${sample}_big_table_filter_${truth_name}_noerror.tsv"

    num_called_noerr=$(grep -E 'UNK|SOM' "${call_peaks_noerr}" | wc -l)
    tp_noerr=$(grep -E 'UNK|SOM' "${both_peaks_noerr}" | awk '{print $4}' | sort | uniq | wc -l) 

    fp_noerr=$((${num_called_noerr} - ${tp_noerr}))
    fn_noerr=$((${num_truth} - ${tp_noerr}))

    precision_noerr=$(echo "scale=4; ${tp_noerr} / (${tp_noerr} + ${fp_noerr})" | bc) 
    recall_noerr=$(echo "scale=4; ${tp_noerr} / (${tp_noerr} + ${fn_noerr})" | bc) 
    f1_noerr=$(echo "scale=4; 2 * (${precision_noerr} * ${recall_noerr}) / (${precision_noerr} + ${recall_noerr})" | bc)

    echo -e "\nAFTER ERROR-PRONE FILTER: ${sample}" >> "${outfile}"
    echo "Number of UNK/SOM peaks called by HAT-seq: ${num_called_noerr}" >> "${outfile}"
    echo "Number of true peaks found by HAT-seq: ${tp_noerr}" >> "${outfile}"
    echo "Precision: ${precision_noerr}" >> "${outfile}"
    echo "Recall: ${recall_noerr}" >> "${outfile}"
    echo "F1-score: ${f1_noerr}" >> "${outfile}"
done
