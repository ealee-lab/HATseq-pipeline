#!/bin/bash
# Sample batchscript to run a simple job on HPC
#SBATCH --time=0-1:00                         # Running time (in hours-minutes-seconds)
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=$1"_HATseq"                         # Job name
#SBATCH --mail-type=FAIL              # send and email when the job begins, ends or fails
#SBATCH --mail-user=shayna.mallett@childrens.harvard.edu          # Email address to send the job status
#SBATCH --nodes=1                               # Number of gpu nodes
#SBATCH --error=./Hatseq_%j.err
#SBATCH --output=./Hatseq_%j.out
#SBATCH -p bch-compute-pe


ls /lab-share/Gene-Lee-ANR-e2 #Rob says this will mount the ANR disk space
#Activate environment 
source activate HATseq_star 
export PATH=~/tools/seqtk:~/tools:~/tools/fastp:$PATH

sample_name=$1
HATseq_version=$2 

snakefile=/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq-pipeline/workflow/Snakefile.${HATseq_version} #_test_rescue #_change_order #_no_blast #main #.test #main
profile=/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq-pipeline/workflow/slurm_carlos
configfile=/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq-pipeline/config/config.*${sample_name}*  #/n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1
directory=/lab-share/Gene-Lee-ANR-e2/shayna/data/output/$1 #/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/$1
#--rerun-triggers mtime
#snakemake  -n  -s ${snakefile}  --configfile ${configfile}  --directory ${directory}
snakemake  --unlock  -s ${snakefile}  --profile ${profile}  --configfile ${configfile}  --directory ${directory}
snakemake   -s ${snakefile}  --profile ${profile}  --configfile ${configfile}  --directory ${directory} 
#snakemake  -s ${snakefile}   -c 4  --configfile ${configfile}  --directory ${directory}
