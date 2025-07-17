#!/bin/bash

# Compare captured peaks to truth set and filter out peaks in error-prone regions.

pdir=$1 # Parent directory containing all pipeline output
error_prone=$2 # Error-prone regions to remove
truth_name=$3 # Name of ground truth set
truth_set=$4 # Path to ground truth set
genome=$5 # Chrom sizes

for call_set in $(ls ${pdir}/*/analysis/*_big_table_filter_reasons.tsv); do
    sample=$(basename "$call_set" _big_table_filter_reasons.tsv)
    outdir=$(dirname "$call_set")
    echo $sample

    truth_out="${outdir}/${sample}_big_table_filter_${truth_name}.tsv"
    error_out="${outdir}/${sample}_big_table_filter_${truth_name}_noerror.tsv"
    class_out="${outdir}/${sample}_big_table_filter_${truth_name}_noerror_unk_som.tsv"

    # Extract relevant fields from peak set then intersect with truth set
    # Fields extracted are chrm, start, end, peak, RPM, strand, classification, filter_reason
    awk -v FS='\t' -v OFS='\t' '{if (NR > 1) print $2,$3,$4,$5,$27,$7,$51,$52}' "$call_set" | \
        bedtools intersect -a stdin \
	    -b <(bedtools slop -b 50 -i "${truth_set}" -g "${genome}") \
	    > "${truth_out}"

    # Remove peaks overlapping error-prone regions
    bedtools intersect -v -a "${truth_out}" \
        -b "${error_prone}" \
	    > "${error_out}"

    # Output file with only UNK and SOM peaks
    grep -E 'UNK|SOM' "${error_out}" > "${class_out}"
done
