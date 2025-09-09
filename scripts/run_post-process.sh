#!/bin/bash

# Post-process each sample by calculating statistics. If your resdir contains
# samples that were not intersected with a truth set, this script will error
# on those samples.
# 
# Usage: bash run_post-process.sh <resdir> <tdir>

#SBATCH --time=0-0:02  # running time (in hours-minutes-seconds)
#SBATCH --cpus-per-task=1  # number of cpus
#SBATCH --mem-per-cpu=80M  # amount of memory (RAM) per cpu
#SBATCH -p bch-compute

resdir=$(realpath $1) # Parent directory with pipeline results
tdir=$(realpath $2) # Directory with truth sets

all_true="${tdir}/l1_hapmapmixture_final_v2_all_tier1-2.bed"
som_true="${tdir}/l1_hapmapmixture_final_v2_somatic_tier1-2.bed"
target_true="${tdir}/l1_hapmapmixture_final_v2_target_tier1-2.bed"
som_target_true="${tdir}/l1_hapmapmixture_final_v2_target_somatic_tier1-2.bed"

for peaks in $(find ${resdir} -name *_candidate_peaks.tsv); do
    sample=$(basename ${peaks} _candidate_peaks.tsv)
    echo $sample
    
    outdir=$(dirname ${peaks})/stats
    mkdir -p ${outdir}
    
    # Calculate precision/recall for each GTS  
    bash calculate_stats.sh \
        "${peaks}" "${all_true}" "${outdir}/${sample}_f1-score_all_tier1-2.txt"

    bash calculate_stats.sh \
        "${peaks}" "${som_true}" "${outdir}/${sample}_f1-score_somatic_tier1-2.txt"

    bash calculate_stats.sh \
        "${peaks}" "${target_true}" "${outdir}/${sample}_f1-score_target_tier1-2.txt"

    bash calculate_stats.sh \
        "${peaks}" "${som_target_true}" "${outdir}/${sample}_f1-score_target_somatic_tier1-2.txt"
done  
