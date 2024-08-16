#!/bin/bash
# Sample batchscript to run a simple job on HPC
#SBATCH --time=3-00:00                         # Running time (in hours-minutes-seconds)
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=run_HATseq                         # Job name
#SBATCH --mail-type=FAIL              # send and email when the job begins, ends or fails
#SBATCH --mail-user=shayna.mallett@childrens.harvard.edu          # Email address to send the job status
#SBATCH --nodes=1                               # Number of gpu nodes
#SBATCH --ntasks=1                              # Number of gpu devices on one gpu node
#SBATCH --error=./Hatseq_%j.err
#SBATCH --output=./Hatseq_%j.out


source activate HATseq_star 
export PATH=~/tools/seqtk:~/tools:~/tools/fastp:$PATH

snakefile=/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/Snakefile.hg38 #_test_rescue #_change_order #_no_blast #main #.test #main
profile=/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/slurm
configfile=/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/config/config.*.$1  #/n/scratch3/users/s/sm745/HATseq_all/config/config.single_cell.$1
directory=/lab-share/Gene-Lee-ANR-e2/shayna/data/output/$1 #/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq_v2/workflow/$1



snakemake --unlock -s ${snakefile}  --profile ${profile}  --configfile ${configfile}  --directory ${directory}
snakemake   -s ${snakefile}  --profile ${profile}  --configfile ${configfile}  --directory ${directory} 


#snakemake  -c 1    -s ${snakefile}   --configfile ${configfile}  --directory ${directory}

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
