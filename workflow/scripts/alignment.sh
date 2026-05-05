#!/usr/bin/env bash
 
exec 2> "${snakemake_log[0]}"

# Map mappable parts of R1 to reference genome
bwa mem -t ${snakemake_params[cores]} -R $(echo ${snakemake_params[RG_ID]}) \
	"${snakemake_input[ref_genome]}" \
	"${snakemake_input[L1HS_trim_fromR1_fq1]}" | \
	samtools view -Sb - \
	> "${snakemake_output[bwa_bam]}"

samtools sort -o "${snakemake_output[bwa_sorted_bam]}" "${snakemake_output[bwa_bam]}" 
samtools index "${snakemake_output[bwa_sorted_bam]}" 

bwa mem -t ${snakemake_params[cores]} -R $(echo ${snakemake_params[RG_ID]}) \
	"${snakemake_input[ref_genome]}" \
	"${snakemake_input[L1HS_polyA_trim_fromR1_fq1]}" | \
	samtools view -Sb - \
	> "${snakemake_output[trim_bwa_bam]}" 
	
# Unique read = read with MAPQ >= 20 and edit distance < 4
samtools view -q 20 "${snakemake_output[trim_bwa_bam]}" | \
	awk '{ for (i=1;i<=NF;i++) \
		{ if ($i ~/AS:i:/) \
			{split ($i,a,":")}; \
		if ($i ~/XS:i:/) \
			{split ($i,b,":")} }; \
		if (a[3]>b[3]) print; }' | \
	awk '{ for (i=1;i<=NF;i++) \
		{ if ($i ~/NM:i:/) \
			{split ($i,a,":"); \
		if (a[3]<4) print; } } }' | \
	cat <(samtools view -H "${snakemake_output[trim_bwa_bam]}") - | \
	samtools view -Sb - \
	> "${snakemake_output[trim_uniq_bam]}"

samtools sort -o "${snakemake_output[trim_uniq_sorted_bam]}" "${snakemake_output[trim_uniq_bam]}" 
samtools index "${snakemake_output[trim_uniq_sorted_bam]}" 

# Find regions (peak) with at least 1 unique read 
cat <( genomeCoverageBed -ibam "${snakemake_output[trim_uniq_sorted_bam]}" -bg -strand + | \
	awk '{OFS="\t"; print $1,$2,$3,$4,".","+",$3-$2}' ) \
	<( genomeCoverageBed -ibam "${snakemake_output[trim_uniq_sorted_bam]}" -bg -strand - | \
	awk '{OFS="\t"; print $1,$2,$3,$4,".","-",$3-$2}' ) | \
	sort -k1,1 -k2,2n | \
	mergeBed -d 1 -i stdin -s -c 1,4,6,7,4 -o count,collapse,distinct,collapse,max | \
	awk -v sample="${snakemake_wildcards[sample]}" \
		'{ OFS="\t"; \
		if($6 == "+") \
			{print $1,$2,$3,sample"-plus-peak-"NR,$4";"$8";"$7";"$5,$6} \
		else \
			{print $1,$2,$3,sample"-minus-peak-"NR,$4";"$8";"$7";"$5,$6} }' \
	>> "${snakemake_output[trim_peaks]}" 

# Output all reads in peak regions
samtools view --region-file "${snakemake_output[trim_peaks]}" \
	"${snakemake_output[bwa_sorted_bam]}" | \
	cat <(samtools view -H "${snakemake_output[bwa_bam]}") - | \
	samtools sort -o "${snakemake_output[peak_sorted_bam]}" - 
 
sambamba index "${snakemake_output[peak_sorted_bam]}"
