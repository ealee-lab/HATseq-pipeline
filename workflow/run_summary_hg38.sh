#!/bin/bash
# Sample batchscript to run a simple job on HPC
#SBATCH --partition=short                         # queue to be used
#SBATCH --time=0-11:00                         # Running time (in hours-minutes-seconds)
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=run_HATseq                         # Job name
#SBATCH --mail-type=FAIL              # send and email when the job begins, ends or fails
#SBATCH --mail-user=shayna.mallett@childrens.harvard.edu          # Email address to send the job status
#SBATCH --output=/n/scratch/users/s/sm745/HATseq_all/scratch/logs/output_%j.txt                  # Name of the output file
#SBATCH --nodes=1                               # Number of gpu nodes
#SBATCH --ntasks=1                              # Number of gpu devices on one gpu node



source activate HATseq_star 
export PATH=~/tools/seqtk:~/tools:$PATH

snakefile=/home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.summary_96cell_hg38 #96cell #Snakefile.combined #.test #main
profile=/home/sm745/HATseq_all/HATseq_v2/workflow/slurm
configfile=/home/sm745/HATseq_all/HATseq_v2/config/config.*.$1  #/n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1
directory=/n/scratch/users/s/sm745/HATseq_all/final/hg38_chr #/n/scratch/users/s/sm745/HATseq_all/analysis_v2_1300/$1



#snakemake  -s ${snakefile}  --profile ${profile}  --configfile ${configfile}  --directory ${directory}
#snakemake  -s ${snakefile} -c 1 --configfile ${configfile}  --directory ${directory} 
snakemake -c 1    -s ${snakefile}   --configfile ${configfile}  --directory ${directory}

#snakemake -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.main --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1



#snakemake --unlock -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module3 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1

#snakemake -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module3 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1



#snakemake --unlock -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module4 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1

#snakemake -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module4 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1

#snakemake --unlock -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module5 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1

#snakemake -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module5 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1

#snakemake --unlock -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module6 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1

#snakemake -s /home/sm745/HATseq_all/HATseq_v2/workflow/Snakefile.module6 --profile /home/sm745/HATseq_all/HATseq_v2/workflow/slurm --configfile /n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1 --directory /n/scratch3/users/s/sm745/HATseq_all/analysis_v2_1300/$1





#snakemake -s /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/Snakefile.module6 --profile /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/slurm --configfile /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/config/config.single_cell --directory /lab-share/Gene-Lee-e2/Public/home/shayna/5666_BA9_1

#snakemake -s /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/Snakefile.module3 --profile /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/slurm --configfile /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/config/config --directory /lab-share/Gene-Lee-e2/Public/home/shayna/5666_BA9_1 

#snakemake -s /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/Snakefile.module2 --profile /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/slurm --configfile /lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/config/config --directory /lab-share/Gene-Lee-e2/Public/home/shayna/5666_BA9_1
