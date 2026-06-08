rule qc:
	input:
		raw_fq1=lambda wc: samples.loc[wc.sample]["raw_fastq1"],
		raw_fq2=lambda wc: samples.loc[wc.sample]["raw_fastq2"],
	output:
		illumina_adapter_trim_fq1=temp("{sample}/{sample}_illumina_adapter_trimmed_R1.fq.gz"),
		illumina_adapter_trim_fq2=temp("{sample}/{sample}_illumina_adapter_trimmed_R2.fq.gz"),
		fastp_html_report="{sample}/reports/{sample}_fastp.html",
		fastp_json_report="{sample}/reports/{sample}_fastp.json",
	params:
		illumina_adapter=config["illumina_adapter"],
		len_req=50,
		report_path=lambda wc: f"{wc.sample}/reports",
	threads: 4
	resources:
		runtime=get_qc_runtime,
		mem_mb=get_qc_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/qc/{sample}.log",
	shell:
		"mkdir -p {params.report_path} \n"
		"fastqc -t {threads} --noextract --nogroup "
			"-o {params.report_path} "
			"{input.raw_fq1} {input.raw_fq2} " # Run fastqc
			"2> {log} \n"

		"fastp -w {threads} -a {params.illumina_adapter} --length_required={params.len_req} "
			"-i {input.raw_fq1} -I {input.raw_fq2} "
			"-o {output.illumina_adapter_trim_fq1} -O {output.illumina_adapter_trim_fq2} "
			"-h {output.fastp_html_report} "
			"-j {output.fastp_json_report} " # Run fastp
			"2>> {log}"

rule preprocessing:
	input:
		illumina_adapter_trim_fq1="{sample}/{sample}_illumina_adapter_trimmed_R1.fq.gz",
		illumina_adapter_trim_fq2="{sample}/{sample}_illumina_adapter_trimmed_R2.fq.gz",
	output:
		L1HS_primer_trim_fromR2_fq1=temp("{sample}/{sample}_L1HS_primer_trimmed_fromR2_R1.fq.gz"),
		L1HS_primer_trim_fromR2_fq2="{sample}/{sample}_L1HS_primer_trimmed_fromR2_R2.fq.gz",
		L1_young_seq_trim_fromR2_fq2=temp("{sample}/{sample}_L1_young_seq_trimmed_fromR2_R2.fq.gz"),
		young_L1_readID_list="{sample}/{sample}_young_L1_readID_list.txt",
		young_L1_reads_fromR1_fq1=temp("{sample}/{sample}_young_L1_reads_fromR1_R1.fq.gz"),
		L1HS_trim_fromR1_fq1="{sample}/{sample}_mappable_L1HS_seq_trimmed_fromR1_R1.fq.gz",
		L1HS_polyA_trim_fromR1_fq1="{sample}/{sample}_mappable_polyA_trimmed_fromR1_R1.fq.gz",
		cutadapt_report="{sample}/reports/{sample}_cutadapt_report.txt",
	params:
		min_len=20,
		L1HS_primer=lambda wc: get_index_seq(wc) + "GGGAGATATACCTAATGCTAGATGACAC",
		L1PA_young_seq="GTTAGTGGGTGCAGCGCACCAGCATGGCACATGTATACATATGTAACTAACCTGCACAATGTGCACATGTACCCTAAAACTTAGAGT",
		L1HSseq="ATTATACTCTAAGTTTTAGGGTACATGTGCACATTGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACCCACTAATGTGTCATCTAGCATTAGGTATATCTCCC",
		adapter_T="TTTTTTTT",
	threads: 4
	resources:
		runtime=get_preprocessing_runtime,
		mem_mb=get_preprocessing_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/preprocessing/{sample}.log",
	script:
		"../scripts/preprocessing.sh"
