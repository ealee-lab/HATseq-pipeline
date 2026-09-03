#!/usr/bin/env bash
 
exec 2> "${snakemake_log[0]}"

# Map mappable parts of R1 to reference genome
bwa mem -t ${snakemake[threads]} -R $(echo ${snakemake_params[RG_ID]}) \
	"${snakemake_input[ref_genome]}" \
	"${snakemake_input[L1HS_trim_fromR1_fq1]}" | \
	samtools view -b - \
	> "${snakemake_output[bwa_bam]}"

samtools sort -o "${snakemake_output[bwa_sorted_bam]}" "${snakemake_output[bwa_bam]}" 
samtools index "${snakemake_output[bwa_sorted_bam]}" 

bwa mem -t ${snakemake[threads]} -R $(echo ${snakemake_params[RG_ID]}) \
	"${snakemake_input[ref_genome]}" \
	"${snakemake_input[L1HS_polyA_trim_fromR1_fq1]}" | \
	samtools view -b - \
	> "${snakemake_output[trim_bwa_bam]}" 
	
# Unique read = read with MAPQ >= 20, AS > XS, and edit distance < 4
samtools view -b -q 20 -e '[AS] > [XS] && [NM] < 4' "${snakemake_output[trim_bwa_bam]}" \
	> "${snakemake_output[trim_uniq_bam]}"

samtools sort -o "${snakemake_output[trim_uniq_sorted_bam]}" "${snakemake_output[trim_uniq_bam]}" 
samtools index "${snakemake_output[trim_uniq_sorted_bam]}" 

# Find regions (peaks) with at least 1 unique read 
cat <( genomeCoverageBed -ibam "${snakemake_input[trim_uniq_sorted_bam]}" -bg -strand + | \
	awk '{OFS="\t"; print $1,$2,$3,$4,".","+",$3-$2}' ) \
	<( genomeCoverageBed -ibam "${snakemake_input[trim_uniq_sorted_bam]}" -bg -strand - | \
	awk '{OFS="\t"; print $1,$2,$3,$4,".","-",$3-$2}' ) | \
	sort -k1,1 -k2,2n | \
	mergeBed -d 1 -i stdin -s -c 1,4,6,7,4 -o count,collapse,distinct,collapse,max | \
	sort -k1,1V -k2,2n | \
	awk -v sample="${snakemake_wildcards[sample]}" \
		'{OFS="\t"; print $1,$2,$3,sample"-uniq-peak-"NR,$4";"$8";"$7";"$5,$6}' \
	>> "${snakemake_output[trim_peaks]}" 

# # Output all reads in peak regions
cat <(samtools view -h -F 16 -L <(awk '{if($6=="+") print}' "${snakemake_output[trim_peaks]}") "${snakemake_input[bwa_sorted_bam]}") \
	<(samtools view -f 16 -L <(awk '{if($6=="-") print}' "${snakemake_output[trim_peaks]}") "${snakemake_input[bwa_sorted_bam]}") | 
	samtools view -b - \
	> "${snakemake_output[peak_bam]}"
 
samtools sort -o "${snakemake_output[peak_sorted_bam]}" "${snakemake_output[peak_bam]}"
samtools index "${snakemake_output[peak_sorted_bam]}"
