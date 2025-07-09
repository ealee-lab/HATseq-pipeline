#!/bin/bash
# Sample batchscript to run a simple job on HPC
#SBATCH --time=0-12:00                         # Running time (in hours-minutes-seconds)
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=HATseq                         # Job name
#SBATCH --mail-type=FAIL              # send and email when the job begins, ends or fails
#SBATCH --mail-user=your.email@you.com          # Email address to send the job status
#SBATCH --nodes=1                               # Number of gpu nodes
#SBATCH --error=./slurm_logs/Hatseq_%j.err
#SBATCH --output=./slurm_logs/Hatseq_%j.out
#SBATCH -p your_partition
#SBATCH -A your_account 


#Activate environment 
source activate HATseq

HATseq_pipeline_directory="/path/to/HATseq-pipeline"
sample_name=$1
HATseq_version=$2 
library=$3
current_dir=$(realpath ./) 
snakefile=${HATseq_pipeline_directory}/workflow/Snakefile.${HATseq_version} 
profile=${HATseq_pipeline_directory}/workflow/slurm
configfile=${HATseq_pipeline_directory}/config/config.${library}.${sample_name}  
directory=/path/to/output/${sample_name} 

snakemake --unlock -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory}
snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} --use-conda --conda-frontend conda
