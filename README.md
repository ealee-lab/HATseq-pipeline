# HATseq-pipeline #
This is a pipeline for the Human Active Trasposon sequencing methodology. 


# Pipeline Setup #
## Conda environment ##
Setting up the environment should require a large amount of memory, best to run on a compute job. Alternatively, you can check the conda.yml file and ensure that the correct versions of all required software are on your PATH and accessible by the pipeline. 
```
cd setup
bash conda_create.sh
conda activate 
```

## Reference Folder ##
```
unzip things, concatenate if needed (or could source original)
```

## Config Files ##
```
cd setup
vim make_config.sh #Edit the paths to the correct files and directories on your server
bash make_config.sh <sample_name> <fastq1> <fastq2> <bulk|microbulk|single> <4_base_in_line_index> [sequencing_platform] [sequencing_center] 
```

# Running the Pipeline # 
```
cd workflow
vim run_HATseq.sh #Edit: your.email@you.com, your_partition, your_account, /path/to/HATseq-pipeline, /path/to/output/
cd slurm
vim cluster-config.yml #Edit: /path/to/HATseq-pipeline, your_partition
bash run_HATseq-pipeline.sh <sample_name> HATseq_hg38 [bulk|micro|single]#Or sbatch to submit as a job on a Slurm HPCC
```

# Demo Data # 
