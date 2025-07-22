#!/bin/bash

# Compare captured peaks to truth set and filter out peaks in error-prone regions.

pdir=$1 # Parent directory containing all pipeline output
error_prone=$2 # Error-prone regions to remove
truth_name=$3 # Name of ground truth set
truth_set=$4 # Path to ground truth set
genome=$5 # Chrom sizes

for call_set in $(ls ${pdir}/*/analysis/*_big_table_filter_reasons.tsv); do
    sample=$(basename "${call_set}" _big_table_filter_reasons.tsv)
    outdir=$(dirname "${call_set}")
    echo $sample

    mkdir -p "${outdir}/post"
    noerr_out="${outdir}/post/${sample}_big_table_filter_noerror.tsv"
    truth_out="${outdir}/post/${sample}_big_table_filter_${truth_name}.tsv"
    truth_noerr_out="${outdir}/post/${sample}_big_table_filter_${truth_name}_noerror.tsv"

    # Extract relevant fields from peak set then intersect
    # Fields extracted are chrm, start, end, peak, RPM, strand, classification, filter_reason
    awk -v FS='\t' -v OFS='\t' '{if (NR > 1) print $2,$3,$4,$5,$27,$7,$51,$52}' "${call_set}" | \
        bedtools intersect -v -a stdin \
        -b "${error_prone}" \
	    > "${noerr_out}"

    awk -v FS='\t' -v OFS='\t' '{if (NR > 1) print $2,$3,$4,$5,$27,$7,$51,$52}' "${call_set}" | \
        bedtools intersect -wo -S -a stdin \
	    -b <(bedtools slop -b 50 -i "${truth_set}" -g "${genome}") \
	    > "${truth_out}"

    # Remove peaks overlapping error-prone regions after intersecting GTS
    bedtools intersect -v -a "${truth_out}" \
        -b "${error_prone}" \
	    > "${truth_noerr_out}"
done
