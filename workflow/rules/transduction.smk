ref_dir = config["reference_dir"]


rule map_transduction_fasta:
	input:
		classified_peaks="{sample}/{sample}_classified_peaks.bed",
		ref_genome=config["bwa_ref_genome"],
		L1HS_decoy=f"{ref_dir}/L1HS_decoy_sequence.fasta",
		softclip3p_fasta="{sample}/{sample}_softclip_3p.fa.gz",
		softclip5p_fasta="{sample}/{sample}_softclip_5p.fa.gz",
		hardclip3p_table="{sample}/{sample}_hardclip_3p.tsv",
		hardclip5p_table="{sample}/{sample}_hardclip_5p.tsv",
	output:
		transduction_bed="{sample}/{sample}_custom_transduction.bed",
		transduction_fasta="{sample}/{sample}_custom_transduction.fasta",
		transduction_fasta_index="{sample}/{sample}_custom_transduction.fasta.bwt",
		softclip3p_trimmed_fasta="{sample}/{sample}_softclip_3p_trimmed.fa.gz",
		softclip5p_trimmed_fasta="{sample}/{sample}_softclip_5p_trimmed.fa.gz",
		softclip3p_trimmed_bam="{sample}/{sample}_softclip_3p_trimmed.bam",
		softclip5p_trimmed_bam="{sample}/{sample}_softclip_5p_trimmed.bam",
		softclip3p_trimmed_bed="{sample}/{sample}_softclip_3p_trimmed.bed",
		softclip5p_trimmed_bed="{sample}/{sample}_softclip_5p_trimmed.bed",
		hardclip3p_bed="{sample}/{sample}_hardclip_3p.bed",
		hardclip5p_bed="{sample}/{sample}_hardclip_5p.bed",
		transduction_report="{sample}/reports/{sample}_transduction_report.txt",
	params:
		adapter_T="TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
		min_length=10,
		RG_ID=get_read_group,
	threads: 4
	resources:
		runtime=get_map_transduction_fasta_runtime,
		mem_mb=get_map_transduction_fasta_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/map_transduction_fasta/{sample}.log",
	group:
		"transduction"
	script:
		"../scripts/map_transduction_fasta.sh"


rule transduction_big_table:
	input:
		peak_readID_list="{sample}/{sample}_peak_readID_list.txt",
		hardclip3p_bed="{sample}/{sample}_hardclip_3p.bed",
		hardclip5p_bed="{sample}/{sample}_hardclip_5p.bed",
		softclip3p_trimmed_bed="{sample}/{sample}_softclip_3p_trimmed.bed",
		softclip5p_trimmed_bed="{sample}/{sample}_softclip_5p_trimmed.bed",
		transduction_bed="{sample}/{sample}_custom_transduction.bed",
	output:
		transduction_big_table="{sample}/{sample}_transduction_big_table.tsv",
	resources:
		runtime=get_transduction_big_table_runtime,
		mem_mb=get_transduction_big_table_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/transduction_big_table/{sample}.log",
	group:
		"transduction"
	script:
		"../scripts/transduction_big_table.py"


rule call_transductions:
	input:
		transduction_bed="{sample}/{sample}_custom_transduction.bed",
		transduction_big_table="{sample}/{sample}_transduction_big_table.tsv",
		classified_peaks="{sample}/{sample}_classified_peaks.bed",
	output:
		transduction_plots="{sample}/{sample}_transduction_plots.pdf",
		filtered_transduction="{sample}/{sample}_filtered_transduction.tsv",
	resources:
		runtime=get_min_runtime,
		mem_mb=get_call_transductions_mem_mb,
	conda:
		"../envs/plotR.yml"
	log:
		"logs/call_transductions/{sample}.log",
	group:
		"transduction"
	shell:
		"Rscript {workflow.basedir}/scripts/transduction_calling.R "
			"{input.transduction_big_table} {output.transduction_plots} "
			"{output.filtered_transduction} {input.transduction_bed} {input.classified_peaks} "
			"2> {log}"
