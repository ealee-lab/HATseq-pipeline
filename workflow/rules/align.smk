rule alignment:
	input:
		ref_genome = config["bwa_ref_genome"],
		L1HS_trim_fromR1_fq1 = "{sample}/{sample}_mappable_L1HS_seq_trimmed_fromR1_R1.fq.gz",
		L1HS_polyA_trim_fromR1_fq1 = "{sample}/{sample}_mappable_polyA_trimmed_fromR1_R1.fq.gz"
	output:
		bwa_bam = temp("{sample}/{sample}_bwa.bam"),
		bwa_sorted_bam = "{sample}/{sample}_bwa_sorted.bam", 
		bwa_sorted_index = "{sample}/{sample}_bwa_sorted.bam.bai",
		trim_bwa_bam = temp("{sample}/{sample}_trim_bwa.bam"),
		trim_uniq_bam = temp("{sample}/{sample}_trim_uniq.bam"), 
		trim_uniq_sorted_bam = temp("{sample}/{sample}_trim_uniq_sorted.bam"),
		trim_uniq_sorted_index = temp("{sample}/{sample}_trim_uniq_sorted.bam.bai"),
		trim_peaks = "{sample}/{sample}_trim_peaks.bed",
		peak_sorted_bam = "{sample}/{sample}_peak_sorted_bwa.bam",
		peak_sorted_index = "{sample}/{sample}_peak_sorted_bwa.bam.bai"
	params:
		RG_ID = get_read_group,
		cores = lambda wc, threads: threads * 2
	threads: 4
	resources:
		runtime = get_alignment_runtime,
		mem_mb = get_alignment_mem_mb
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/alignment/{sample}.log"
	script:
		"../scripts/alignment.sh"
		