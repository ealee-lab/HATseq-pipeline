#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"

intersectBed -wao -a "${snakemake_input[peaks]}" -b "${snakemake_input[satellites]}" | \
	mergeBed -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	> "${snakemake_output[satellite_intersect]}"

intersectBed -wao -a "${snakemake_input[peaks]}" -b "${snakemake_input[segdups]}" | \
	mergeBed -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	> "${snakemake_output[segdup_intersect]}"

intersectBed -wa -wb \
	-a <(awk -v OFS='\t' '{if($6 == "+") print $1,$3-30,$3,$4,$5,$6}' "${snakemake_input[peaks]}") \
	-b <(grep "Homo_T" "${snakemake_input[homopolymers]}") | \
	mergeBed -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	> "${snakemake_output[homopolymer_intersect]}"

intersectBed -wa -wb \
	-a <(awk -v OFS='\t' '{if($6 == "-") print $1,$2,$2+30,$4,$5,$6}' "${snakemake_input[peaks]}") \
	-b <(grep "Homo_A" "${snakemake_input[homopolymers]}") | \
	mergeBed -s -i stdin -c 4,10 -o distinct,collapse | \
	cut -f4,5 \
	>> "${snakemake_output[homopolymer_intersect]}"
