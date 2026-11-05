#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"

bedtools intersect -wao -a "${snakemake_input[peaks]}" -b "${snakemake_input[satellites]}" | \
	bedtools merge -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	> "${snakemake_output[satellite_intersect]}"

bedtools intersect -wao -a "${snakemake_input[peaks]}" -b "${snakemake_input[segdups]}" | \
	bedtools merge -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	> "${snakemake_output[segdup_intersect]}"

bedtools intersect -wa -wb \
	-a <(awk -v OFS='\t' '{if($6 == "+") print}' "${snakemake_input[breakends]}") \
	-b <(grep "Homo_T" "${snakemake_input[homopolymers]}") | \
	bedtools merge -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	> "${snakemake_output[homopolymer_intersect]}"

bedtools intersect -wa -wb \
	-a <(awk -v OFS='\t' '{if($6 == "-") print}' "${snakemake_input[breakends]}") \
	-b <(grep "Homo_A" "${snakemake_input[homopolymers]}") | \
	bedtools merge -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	>> "${snakemake_output[homopolymer_intersect]}"
