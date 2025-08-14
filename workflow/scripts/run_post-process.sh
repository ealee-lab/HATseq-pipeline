#!/bin/bash

# Post-process each sample by calculating statistics.

resdir=$(realpath $1) # Parent directory with pipeline results
tdir=$(realpath $2) # Directory with truth sets

all_true="${tdir}/l1_hapmapmixture_final_v2_all_tier1-2.bed"
som_true="${tdir}/l1_hapmapmixture_final_v2_somatic_tier1-2.bed"
target_true="${tdir}/l1_hapmapmixture_final_v2_target_tier1-2.bed"
som_target_true="${tdir}/l1_hapmapmixture_final_v2_target_somatic_tier1-2.bed"

for path in $(ls ${resdir}/); do
    sample=$(basename "${path}")
    echo $sample
    
    indir="${resdir}/${sample}/analysis"
    outdir="${resdir}/${sample}/analysis/stats"
    mkdir -p "${outdir}"
    
    big_table="${indir}/${sample}_filtered_peaks.tsv"

    # Calculate precision/recall for each GTS  
    bash calculate_stats.sh \
        "${big_table}" "${all_true}" "${outdir}/${sample}_f1-score_all_tier1-2.txt"

    bash calculate_stats.sh \
        "${big_table}" "${som_true}" "${outdir}/${sample}_f1-score_somatic_tier1-2.txt"

    bash calculate_stats.sh \
        "${big_table}" "${target_true}" "${outdir}/${sample}_f1-score_target_tier1-2.txt"

    bash calculate_stats.sh \
        "${big_table}" "${som_target_true}" "${outdir}/${sample}_f1-score_target_somatic_tier1-2.txt"
done  
