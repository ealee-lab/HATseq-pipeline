# HAT-seq Pipeline #
This is a pipeline for the Human Active Transposon sequencing methodology. 

# Setting up the Pipeline #
## Reference Folder ##
This pipeline uses several annotations to remove reference and non-reference L1 elements from the set of putative somatic insertions. Error-prone regions are also used to identify regions prone to noise. This folder contains the necessary annotation files and error-prone regions as listed:

```
- reference
    - 1kgp
        - ALL_MELT_ME_1000G_HC_20190901.AF.bed # 1000G_2504_high_coverage_SV collection from 1000 Genomes
    - 1019_ONT
        - 1019_ONT_Schloissnig2025.bed # Schloissnig et al. (2025), Supplementary Table 18
    - ErrorProne
        - bulk.NIH_Aging-AT.merged.error-prone.bed # generated from NIH NeuroBioBank data
        - bulk.NIH_Aging.merged.error-prone.bed # generated from NIH NeuroBioBank data
    - gnomAD-SV
        - gnomad.v4.1.ME.sites.bed # gnomAD v4.1
    - HGSVC3
        - MEI_Callset_GRCh38.ALL.20241211.bed # Logsdon et al. (2025), Supplementary Table 32
        - Ortho_MEI_GRCh38.ALL.20241211.bed # Logsdon et al. (2025), Supplementary Table 34
    - human
        - hg38.genome # list of chrom sizes
        - hg38.hg19.homo8.chr.bed # homopolymers
    - nyuwa
        - MEI.GRCh38.HMEIDv1.1.final.bed # NyuWa
    - RepeatMasker
        - hg38.Evrony_KR_960.v37.revised.bed # Evrony et al. (2012)
        - hg38.repeatmasker.bed # RepeatMasker 
    - SegDup
        - hg38.genomicSuperDups.v37.chr.bed # segmental duplications
    - xTea
        - xTea_Borges-Monroy2021.accessioned.bed # Chu et al. (2021)
    - L1HS_decoy_sequence.fasta
```

## Config Files ##
See `config.library.yml` for an example config file. It is expected that you will have one config file for each library type (bulk, microbulk, single-cell), but you can split the pipeline into as many runs as you would like. Each config file links to a sample sheet with sample-specific information (`samples.library.tsv`) including donor, tissue/cell/rep, raw_fastq1, raw_fastq2, and index_seq. Users can compare across either tissues, cells, replicates, or both tissues and replicates. Replicates must be identified by a capital letter followed by a number (e.g. R1) in the sample sheet.

The error_prone and truth_set fields in the config file are optional. You can use the provided error-prone regions or generate your own with the `run_make_error-prone_regions.sh` script. The truth_set field is available in the case that you want to compare the pipeline output to a ground truth set of insertions. Truth sets are expected to be BED6 files of L1 insertions.

```
cd config
cp samples.library.tsv <your_sample_sheet>
vim <your_sample_sheet> # add your sample information

cp config.library.yml <your_config_file>
vim <your_config_file> # edit the parameters 
```

# Running the Pipeline # 
You must have conda installed to run this pipeline as it will install conda environments as needed. Also, this pipeline uses a profile to run as a job on a Slurm HPCC. If running locally, you can specify your own custom configuration profile. 

```
cd workflow/slurm
vim cluster-config.yml # edit /path/to/HATseq-pipeline, your_partition

cd .. # should be in workflow directory
vim run_HATseq-pipeline.sh # edit Slurm parameters (i.e. the SBATCH lines)
sbatch run_HATseq-pipeline.sh <config> <outdir> <tmpdir> <profile> # submit job on a Slurm HPCC
```

Note that this pipeline was tested with Snakemake 7.32.4, so changes may be necessary to make it compatible with newer Snakemake versions.

# Debugging the Pipeline # 
One of the most common issues you may come across while running the pipeline are Out of Memory or Timeout errors. Each step in the pipeline has set resource limits based on testing. These limits are implemented so as to optimize job efficiency on HPC clusters, especially when running many samples. 

That being said, if the allocated resources are insufficient, users have several options. Increasing `restart-times` in `workflow/slurm/config.yaml` will increase the number of times rules can be restarted, with each attempt being allocated double the resources of the previous attempt. 

For a more permanent solution, users should locate the corresponding rule and resource in `rules/common.smk` and increase resources. For example, if the pipeline errors on the preprocessing rule with an Out of Memory error, locate the `get_preprocessing_mem_mb()` function in `common.smk` and type in amount of memory (in MB) you would like to request for that rule.

Lastly, users can override the resource requirements set in `common.smk` by adding the options `--default-resources RESOURCE=VALUE` or `--set-resources RULE:RESOURCE=VALUE` (for specific rules) to the running command in `run_HATseq-pipeline.sh`.
