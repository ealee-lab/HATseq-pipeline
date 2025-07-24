#!/bin/bash

# Post-process each sample by filtering out error-prone regions, comparing to truth set, 
# and calculating statistics.

resdir=$(realpath $1) # Parent directory with pipeline results
tdir=$(realpath $2) # Directory with truth sets
err=$3 # Error-prone regions to remove
genome=$4 # Chrom sizes

all_true="${tdir}/l1_hapmapmixture_final_v2_all_tier1-2.bed"
som_true="${tdir}/l1_hapmapmixture_final_v2_somatic_tier1-2.bed"
target_true="${tdir}/l1_hapmapmixture_final_v2_target_tier1-2.bed"
som_target_true="${tdir}/l1_hapmapmixture_final_v2_target_somatic_tier1-2.bed"

# Find lists of insertions to remove for each subset GTS
# e.g. If comparing to somatic GTS, remove germline GTS insertions
grep germline "${all_true}" | \
    awk -v FS='\t' '{print $4}' \
    > "${tdir}/l1_hapmapmixture_final_v2_germline_ins.txt"

awk -v FS='\t' '{print $4}' "${target_true}" | \
    grep -v -f - "${all_true}" | \
    awk -v FS='\t' '{print $4}' \
    > "${tdir}/l1_hapmapmixture_final_v2_non-target_ins.txt"

cat "${tdir}/l1_hapmapmixture_final_v2_germline_ins.txt" "${tdir}/l1_hapmapmixture_final_v2_non-target_ins.txt" \
    > "${tdir}/l1_hapmapmixture_final_v2_germline_non-target_ins.txt"

for path in $(ls ${resdir}/); do
    sample=$(basename "${path}")
    echo $sample
    
    indir="${resdir}/${sample}/analysis"
    outdir="${resdir}/${sample}/analysis/post"
    mkdir -p "${outdir}"
    
    big_table="${indir}/${sample}_big_table_filter_reasons.tsv"
    out_file="${outdir}/${sample}_big_table_filter_all_tier1-2.tsv"
    noerr_file="${outdir}/${sample}_big_table_filter_all_tier1-2_noerror.tsv"

    # Compare HAT-seq peaks to each GTS and filter out peaks in error-prone regions
    bash compare_truth_and_filter_error.sh \
        "${big_table}" "${out_file}" "${noerr_file}" "${all_true}" "${err}" "${genome}"

    # Calculate precision/recall for each GTS  
    bash calculate_stats.sh \
        "${out_file}" "${noerr_file}" "${outdir}/${sample}_f1-score_all_tier1-2.txt" "${all_true}" 

    bash calculate_stats.sh \
        "${out_file}" "${noerr_file}" "${outdir}/${sample}_f1-score_somatic_tier1-2.txt" "${som_true}" \
        "${tdir}/l1_hapmapmixture_final_v2_germline_ins.txt"

    bash calculate_stats.sh \
        "${out_file}" "${noerr_file}" "${outdir}/${sample}_f1-score_target_tier1-2.txt" "${target_true}" \
        "${tdir}/l1_hapmapmixture_final_v2_non-target_ins.txt"

    bash calculate_stats.sh \
        "${out_file}" "${noerr_file}" "${outdir}/${sample}_f1-score_target_somatic_tier1-2.txt" "${som_target_true}" \
        "${tdir}/l1_hapmapmixture_final_v2_germline_non-target_ins.txt"
done  
