#!/usr/bin/env bash

set -x
exec 2> "${snakemake_log[0]}"

# Remove L1HS adapter from R2, discarding pairs where R2 does not contain adapter
# R1 and R2 fastqs are swapped and the sequence is with respect to the 5' of R2
cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len]} \
	-g ${snakemake_params[adapter_L1HS]} \
	--discard-untrimmed -O 10 -e 0.08 \
	-p "${snakemake_output[L1HS_adapter_trimmed_fromR2_fastq1]}" \
	-o "${snakemake_output[L1HS_adapter_trimmed_fromR2_fastq2]}" \
	"${snakemake_input[illumina_adapter_trimmed_fastq2]}" \
	"${snakemake_input[illumina_adapter_trimmed_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"

# Remove L1HS primer from R2, discarding pairs where R2 does not contain primer
# R1 and R2 fastqs are swapped and the sequence is with respect to the 5' of R2
cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len]} \
	-g ${snakemake_params[adapter_primer]} \
	--discard-untrimmed -O 5 -e 0.16 \
	-p "${snakemake_output[L1HS_primer_trimmed_fromR2_fastq1]}" \
	-o "${snakemake_output[L1HS_primer_trimmed_fromR2_fastq2]}" \
	"${snakemake_output[L1HS_adapter_trimmed_fromR2_fastq2]}" \
	"${snakemake_output[L1HS_adapter_trimmed_fromR2_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"  

# Remove R2s that don't contain L1PA-specific sequence
cutadapt -j ${snakemake[threads]} \
	-g ${snakemake_params[L1PA_young_seq]} \
	--discard-untrimmed -e 0.05 -O 20 \
	"${snakemake_output[L1HS_adapter_trimmed_fromR2_fastq2]}" | \
	seqkit seq -i -n - \
	> "${snakemake_output[subfamily_specific_R2_list]}"

# Remove the L1HS adapter from R1
cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len]} \
	-a ${snakemake_params[adapter_overlap_L1HS]} \
	-O 5 -e 0.08 \
	-o "${snakemake_output[L1HS_adapter_trimmed_fromR1_fastq1]}" \
	"${snakemake_output[L1HS_primer_trimmed_fromR2_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"

# Remove the L1HS primer from R1
cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len]} \
	-a ${snakemake_params[adapter_overlap_primer]} \
	-O 5 -e 0.16 \
	--untrimmed-output "${snakemake_output[L1HS_untrimmed_fromR1_fastq1]}" \
	-o "${snakemake_output[L1HS_primer_trimmed_fromR1_fastq1]}" \
	"${snakemake_output[L1HS_adapter_trimmed_fromR1_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"  

cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len_informative]} \
	"${snakemake_output[L1HS_primer_trimmed_fromR1_fastq1]}" | \
	seqtk subseq - "${snakemake_output[subfamily_specific_R2_list]}" | \
	seqkit seq -i -n -  \
	> "${snakemake_output[mappable_R1_list]}" 

# Create list of R1s that have a pair that passed 
seqtk subseq "${snakemake_output[L1HS_untrimmed_fromR1_fastq1]}" \
	"${snakemake_output[subfamily_specific_R2_list]}" | \
	seqkit seq -i -n -  \
	>> "${snakemake_output[mappable_R1_list]}"  

seqtk subseq "${snakemake_output[L1HS_adapter_trimmed_fromR2_fastq1]}" \
	"${snakemake_output[mappable_R1_list]}" | \
	gzip - \
	> "${snakemake_output[mappable_R1_fastq1]}"  

cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len]} \
	-a ${snakemake_params[L1HSseq]} \
	-O ${snakemake_params[min_ovlp]} -e ${snakemake_params[max_err]} \
	-o "${snakemake_output[mappable_R1_L1HSseq_trimmed_fromR1_fastq1]}" \
	"${snakemake_output[mappable_R1_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}" 

cutadapt -j ${snakemake[threads]} \
	-a ${snakemake_params[adapter_T]} -e ${snakemake_params[max_err_T]} \
	-o "${snakemake_output[mappable_R1_L1HSseq_polyA_trimmed_fromR1_fastq1]}" \
	"${snakemake_output[mappable_R1_L1HSseq_trimmed_fromR1_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"
