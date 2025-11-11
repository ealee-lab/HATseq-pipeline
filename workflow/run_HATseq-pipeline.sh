#!/bin/bash

# Use Slurm to run the HATseq pipeline

#SBATCH --time=0-12:00  # running time (in hours-minutes-seconds)
#SBATCH --cpus-per-task=1  # number of cpus
#SBATCH --mem-per-cpu=1G  # amount of memory (RAM) per cpu
#SBATCH --job-name=HATseq  # job name
#SBATCH --mail-type=FAIL   # send and email when the job fails
#SBATCH --mail-user=your.email@you.com  # email address to send the job status
#SBATCH --nodes=1  # number of gpu nodes
#SBATCH --error=./slurm_logs/Hatseq_%j.err
#SBATCH --output=./slurm_logs/Hatseq_%j.out
#SBATCH --error=./slurm_logs/Hatseq_%j.err
#SBATCH -p bch-compute
 
config=$(realpath $1)
outdir=$(realpath $2)
tmpdir=$(realpath $3)
profile=$(realpath $4)

# Create environment if needed and activate 
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

if conda env list | grep -E "^HATseq\b"; then
    conda activate HATseq
else
    conda env create --file envs/conda.yml --yes
    conda activate HATseq
fi

export TMPDIR=${tmpdir} # used as tmpdir by snakemake

# Run the pipeline
snakemake --unlock -s Snakefile \
    --configfile ${config} \
    --directory ${outdir} \
    --profile ${profile}

snakemake -s Snakefile \
    --configfile ${config} \
    --directory ${outdir} \
    --profile ${profile} \
    --use-conda --conda-frontend conda
