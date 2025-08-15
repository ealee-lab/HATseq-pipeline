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
#SBATCH -p your_partition
#SBATCH -A your_account 

snakefile=$1 
profile=$2
config=$3
outdir=$4
tmpdir=$5

# Create enviornment if needed and activate 
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

if conda env list | grep -E "^HATseq\b"; then
    echo "here"
    conda activate HATseq
else
    conda env create --file envs/conda.yml --yes
    conda activate HATseq
fi

export TMPDIR=${tmpdir} # used as tmpdir by snakemake

# Run the pipeline
snakemake --unlock -s ${snakefile} \
    --profile ${profile} \
    --configfile ${config} \
    --directory ${outdir}

snakemake -s ${snakefile} \
    --profile ${profile} \
    --configfile ${config} \
    --directory ${outdir} \
    --use-conda --conda-frontend conda
