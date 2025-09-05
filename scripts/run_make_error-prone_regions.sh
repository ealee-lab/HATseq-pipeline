#!/bin/bash
# Make error-prone regions from recurrent noise peaks

#SBATCH --time=0-0:01  # running time (in hours-minutes-seconds)
#SBATCH --cpus-per-task=1  # number of cpus
#SBATCH --mem-per-cpu=80M  # amount of memory (RAM) per cpu
#SBATCH --job-name=errorprone  # job name
##SBATCH --mail-type=END  # send and email when the job ends
##SBATCH --mail-user=rohini.gadde@childrens.harvard.edu  # email address to send the job status
#SBATCH --error=./slurm_logs/errorprone_%j.err
#SBATCH --output=./slurm_logs/errorprone_%j.out
#SBATCH -p bch-compute

donor_list=$(realpath $1)
resdir=$(realpath $2)
outdir=$(realpath $3)

# Activate environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate HATseq

python make_error-prone_regions.py ${donor_list} ${resdir} ${outdir}
