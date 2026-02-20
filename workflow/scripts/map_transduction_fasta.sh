#!/usr/bin/env bash

set -x 
exec 2> "${snakemake_log[0]}"

awk -v OFS='\t' 'NR>1 {if ($7 != "FP") print $1,$2,$3,$4}' "${snakemake_input[classified_peaks]}" | \
	sort -k1,1 -k2,2n | \
	mergeBed -i stdin -c 4 -o collapse | \
	awk -v sample="${snakemake_wildcards[sample]}" -v OFS='\t' \
		'{print $1,$2,$3,sample"_custom_reference_"NR,$4}' \
	> "${snakemake_output[transduction_bed]}"

# Everything is with respect to the reference at this point
bedtools getfasta -fi "${snakemake_input[ref_genome]}" \
	-bed "${snakemake_output[transduction_bed]}" -name \
	> "${snakemake_output[transduction_fasta]}"  

cat "${snakemake_input[L1HS_decoy]}" \
	>> "${snakemake_output[transduction_fasta]}"

bwa index "${snakemake_output[transduction_fasta]}" 

cutadapt -g X${snakemake_params[adapter_T]} -m ${snakemake_params[min_length]} \
	-o "${snakemake_output[softclip3p_trimmed_fasta]}" \
	"${snakemake_input[softclip3p_fasta]}" \
	>> "${snakemake_output[transduction_report]}"

cutadapt -a ${snakemake_params[adapter_T]}X -m ${snakemake_params[min_length]} \
	-o "${snakemake_output[softclip5p_trimmed_fasta]}" \
	"${snakemake_input[softclip5p_fasta]}" \
	>> "${snakemake_output[transduction_report]}" 

bwa mem -t ${snakemake[threads]} -R $(echo ${snakemake_params[RG_ID]}) \
	"${snakemake_output[transduction_fasta]}" "${snakemake_output[softclip3p_trimmed_fasta]}" | \
	samtools sort -o ${snakemake_output[softclip3p_trimmed_bam]}  
bwa mem -t ${snakemake[threads]} -R $(echo ${snakemake_params[RG_ID]}) \
	"${snakemake_output[transduction_fasta]}" "${snakemake_output[softclip5p_trimmed_fasta]}" | \
	samtools sort -o "${snakemake_output[softclip5p_trimmed_bam]}" 

samtools index "${snakemake_output[softclip3p_trimmed_bam]}" 
samtools index "${snakemake_output[softclip5p_trimmed_bam]}" 

samtools view -q 20 "${snakemake_output[softclip3p_trimmed_bam]}" | \
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
	cat <(samtools view -H "${snakemake_output[softclip3p_trimmed_bam]}") - | \
	samtools view -Sb - | \
	bamToBed -i stdin \
	> "${snakemake_output[softclip3p_trimmed_bed]}"
  
samtools view -q 20 "${snakemake_output[softclip5p_trimmed_bam]}" | \
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
	cat <(samtools view -H "${snakemake_output[softclip5p_trimmed_bam]}") - | \
	samtools view -Sb - | \
	bamToBed -i stdin \
	> "${snakemake_output[softclip5p_trimmed_bed]}" 

intersectBed -wa -wb \
	-a <(awk '{OFS="\t"; split($2,a,":"); split(a[2],b,"-"); print a[1],b[1],b[2],$1,".",a[3]}' "${snakemake_input[hardclip3p_table]}") \
	-b "${snakemake_output[transduction_bed]}" \
	> "${snakemake_output[hardclip3p_bed]}"

intersectBed -wa -wb \
	-a <(awk '{OFS="\t"; split($2,a,":"); split(a[2],b,"-"); print a[1],b[1],b[2],$1,".",a[3]}' "${snakemake_input[hardclip5p_table]}") \
	-b "${snakemake_output[transduction_bed]}" \
	> "${snakemake_output[hardclip5p_bed]}" 		
