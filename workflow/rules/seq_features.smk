rule gmotif:
	input:
		L1HS_primer_trim_fromR2_fq2="{sample}/{sample}_L1HS_primer_trimmed_fromR2_R2.fq.gz",
		young_L1_readID_list="{sample}/{sample}_young_L1_readID_list.txt",
	output:
		pass_gmotif_list="{sample}/{sample}_pass_gmotif_list.txt",
	resources:
		runtime=get_gmotif_runtime,
		mem_mb=get_gmotif_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/gmotif/{sample}.log",
	shell:
		"seqtk subseq {input.L1HS_primer_trim_fromR2_fq2} {input.young_L1_readID_list} | "
			"seqtk seq -a - | "
			"grep 'CTTAGAGT' -B 1 | "
			"seqkit seq -i -n - | sort | uniq "
			"> {output.pass_gmotif_list} 2> {log}"


rule extract_clipping:
	input:
		peak_sorted_bam="{sample}/{sample}_peak_sorted_bwa.bam",
	output:
		softclip3p_fasta="{sample}/{sample}_softclip_3p.fa.gz",
		softclip3p_table = "{sample}/{sample}_softclip_3p.tsv",
		softclip5p_fasta="{sample}/{sample}_softclip_5p.fa.gz",
		hardclip3p_table="{sample}/{sample}_hardclip_3p.tsv",
		hardclip5p_table="{sample}/{sample}_hardclip_5p.tsv",
	resources:
		runtime=get_extract_clipping_runtime,
		mem_mb=get_extract_clipping_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/extract_clipping/{sample}.log",
	script:
		"../scripts/extract_clipping.py"


# rule extract_endpoints:
# 	input:
# 		peak_sorted_bam="{sample}/{sample}_peak_sorted_bwa.bam",
# 	output:
# 		endpoint3p_table="{sample}/{sample}_endpoint_3p.tsv",
# 	resources:
# 		runtime=get_extract_clipping_runtime,
# 		mem_mb=get_extract_endpoints_mem_mb,
# 	conda:
# 		"../envs/HATseq.yml"
# 	log:
# 		"logs/extract_endpoints/{sample}.log",
# 	script:
# 		"../scripts/extract_endpoints.py"


rule junction_spanning:
	input:
		softclip3p_table="{sample}/{sample}_softclip_3p.tsv",
		# uniq_readID_list="{sample}/{sample}_unique_readID_list.txt"
	output:
		polyT_reads="{sample}/{sample}_polyT_reads.txt",
	params:
		adapter_T="TTTTTTTT",
	resources:
		runtime=get_min_runtime,
		mem_mb=get_min_mem_mb,
	log:
		"logs/junction_spanning/{sample}.log",
	shell:
		"cat {input.softclip3p_table} | "
			"awk -v OFS='\t' "
			"'{{if(substr($5,1,30) ~ \"{params.adapter_T}\") print $1,$2,$3,$4,substr($5,1,30)}}' "
			"> {output.polyT_reads} 2> {log}"
