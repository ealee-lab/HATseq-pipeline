#!/bin/bash

# Compare captured peaks to truth set and filter out peaks in error-prone regions.

big_table=$1 # Path to big_table_filter_reasons file
out_file=$2 # Path to output file intersected with ground truth set
noerr_file=$3 # Path to output file with error-prone regions removed
truth_set=$4 # Path to ground truth set
error_prone=$5 # Error-prone regions to remove
genome=$6 # Chrom sizes

# Extract relevant fields from peak set then intersect with truth set
# Fields extracted are chrm, start, end, peak, RPM, strand, classification, filter_reason
awk -v FS='\t' -v OFS='\t' '{if (NR > 1) print $2,$3,$4,$5,$27,$7,$51,$52}' "${big_table}" | \
    bedtools intersect -wao -S -a stdin \
    -b <(bedtools slop -b 50 -i "${truth_set}" -g "${genome}") \
    > "${out_file}"

# Remove peaks overlapping error-prone regions after intersecting GTS
bedtools intersect -v -a "${out_file}" \
    -b "${error_prone}" \
    > "${noerr_file}"
