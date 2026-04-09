#!/usr/bin/env bash

set -x
exec 2> "${snakemake_log[0]}"

# Remove L1HS AC primer from R2 
# Discard pairs if R2 does not contain sequence
cutadapt -j ${snakemake[threads]} \
	-g "${snakemake_params[L1HS_primer]}" \
	--discard-untrimmed -O 10 -e 0.08 \
	-o "${snakemake_output[L1HS_primer_trimmed_fromR2_fastq2]}" \
	-p "${snakemake_output[L1HS_primer_trimmed_fromR2_fastq1]}" \
	"${snakemake_input[illumina_adapter_trimmed_fastq2]}" \
	"${snakemake_input[illumina_adapter_trimmed_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"

# Make list of L1HS-derived read IDs
cutadapt -j ${snakemake[threads]} \
	-g "^${snakemake_params[L1PA_young_seq]}" \
	--discard-untrimmed -e 4 \
	-o "${snakemake_output[L1_young_seq_trimmed_fromR2_fastq2]}" \
	"${snakemake_output[L1HS_primer_trimmed_fromR2_fastq2]}" \
	>> "${snakemake_output[cutadapt_report]}" 

seqkit seq -i -n "${snakemake_output[L1_young_seq_trimmed_fromR2_fastq2]}" \
	> "${snakemake_output[young_L1_readID_list]}"

# Output all R1s whose mate is an L1HS-derived read
seqtk subseq "${snakemake_output[L1HS_primer_trimmed_fromR2_fastq1]}" \
	"${snakemake_output[young_L1_readID_list]}" | \
	gzip - \
	> "${snakemake_output[young_L1_reads_fromR1_fastq1]}"

# Trim L1HS sequence from R1 if present
# If present, cutadapt will also remove any tailing primer sequences
cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len]} \
	-a "${snakemake_params[L1HSseq]}" \
	-O 20 -e 0.08 \
	-o "${snakemake_output[mappable_L1HSseq_trimmed_fromR1_fastq1]}" \
	"${snakemake_output[young_L1_reads_fromR1_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"

# Trim polyA tail from R1 if present
cutadapt -j ${snakemake[threads]} \
	-m ${snakemake_params[min_len]} \
	-a "${snakemake_params[adapter_T]}" \
	-O 6 -e 0 \
	-o "${snakemake_output[mappable_L1HSseq_polyA_trimmed_fromR1_fastq1]}" \
	"${snakemake_output[mappable_L1HSseq_trimmed_fromR1_fastq1]}" \
	>> "${snakemake_output[cutadapt_report]}"
