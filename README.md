# HATseq-pipeline #
This is a pipeline for the Human Active Transposon sequencing methodology. 


# Pipeline Setup #
## Conda environment ##
Setting up the HATseq environment should require a larger amount of memory, best to run on a compute job if running on an HPC cluster. Note that you must already have conda installed before this step.

```
cd setup
conda env create --file conda.yml
```

## Reference Folder ##
This folder contains the following annotation files:
```
- reference
    - 1kgp
        - 1kgp.AF.ME.final.hg38.bed
    - gnomAD-SV
        - gnomad_v4_MEI.sites_hg38.revised.chr.bed
    - human
        - hg38.hg19.homo8.chr.bed # homopolymers
        - /human/hg38.genome # list of chrom sizes
    - nyuwa
        - nyuwa/MEI.GRCh38.HMEIDv1.1.final.bed
    - RepeatMasker
        - hg38.Evrony_KR_960.v37.revised.bed.chr.bed
        - repeatmasker.v38.bed
    - SegDup
        - hg38.genomicSuperDups.v37.chr.bed
    - xTea
        - xTea-Borges-Monroy2021.accessioned.bed
    - L1HS_decoy_sequence.fasta
```

## Config Files ##
See `config.library.yml` for an example config file. It is expected that you will have one config file for each library type (bulk, microbulk, single-cell), but you can split the pipeline into as many runs as you would like. Each config file links to a sample sheet with sample-specific information (`samples.library.tsv`) including sample_name, raw_fastq1, raw_fastq2, index_seq, error_prone, and truth_set. Each sample can have different error_prone and truth_set files if desired or the field can be set to "-NA-" if not desired.

```
unzip things, concatenate if needed (or could source original)
```

## Config Files ##
```
cd config
cp samples.library.tsv <your_sample_sheet>
vim <your_sample_sheet> # add your sample information

cp config.library.yml <your_config_file>
vim <your_config_file> # edit the parameters 
```

# Running the Pipeline # 
This pipeline uses a profile to run the pipeline as a job on a Slurm HPCC. 
```
cd workflow/slurm
vim cluster-config.yml # edit /path/to/HATseq-pipeline, your_partition

cd .. # should be in workflow directory
vim run_HATseq-pipeline.sh # edit Slurm parameters (SBATCH lines)
bash run_HATseq-pipeline.sh <snakefile> <your_config_file> <profile> <results_dir> <temp_dir> # or sbatch to submit as a job on a Slurm HPCC
```
