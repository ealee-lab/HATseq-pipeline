rule create_big_table:
	input:
		peaks = "{sample}/{sample}_peaks.bed",
		peak_readID_list = "{sample}/{sample}_peak_readID_list.txt",
		intersect_annotated = "{sample}/{sample}_intersect_annotated.bed",
		pass_gmotif_list = "{sample}/{sample}_pass_gmotif_list.txt",
		# Ntag_list = "{sample}/{sample}_Ntag_list.txt",
		unique_start_positions = "{sample}/{sample}_unique_start_positions.bed",
		nearby_peaks = "{sample}/{sample}_nearby_peaks.tsv",
		satellite_intersect="{sample}/{sample}_satellite_intersect.bed",
        segdup_intersect="{sample}/{sample}_segdup_intersect.bed",
        homopolymer_intersect="{sample}/{sample}_homopolymer_intersect.bed",
		max_depth = "{sample}/{sample}_max_depth_distance_to_boundary.txt",
		polyT_reads = "{sample}/{sample}_polyT_reads.txt",
		endpoint3p_table = "{sample}/{sample}_endpoint_3p.tsv"
	output:
		big_table = "{sample}/{sample}_big_table.tsv"
	resources:
		runtime = get_create_big_table_runtime,
		mem_mb = get_create_big_table_mem_mb
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/create_big_table/{sample}.log"
	group: "classify"
	script:
		"../scripts/create_big_table.py"

rule filter_and_classify:
	input:
		big_table = "{sample}/{sample}_big_table.tsv",
	output:
		filter_reasons = "{sample}/{sample}_big_table_filter_reasons.tsv",
		classified_peaks = "{sample}/{sample}_classified_peaks.bed",
		plots = "{sample}/{sample}_filtering_plots.pdf",
		summary_file = "{sample}/{sample}_summary.txt",
		igv_script = "{sample}/{sample}_igv.sh"
	params:
		library = config["library"],
		error_prone = config.get("error_prone", "-NA-"),
		truth_set = config.get("truth_set", "-NA-"), 
		path = lambda wc: os.getcwd() + f"/{wc.sample}/igv/"
	resources: 
		runtime = get_min_runtime,
		mem_mb = get_filter_and_classify_mem_mb,
	conda:
		"../envs/plotR.yml"
	log:
		"logs/filter_and_classify/{sample}.log"
	group: "classify"
	shell:
		"Rscript {workflow.basedir}/scripts/filter_and_plot.R "
			"{input.big_table} {output.filter_reasons} "
			"{params.library} {output.plots} "
			"{output.classified_peaks} {output.summary_file} "
			"{params.error_prone} {params.truth_set} "
			"2> {log} \n"

		"grep -E 'UNK|SOM' {output.classified_peaks} | "
			"awk -v OFS='\t' '{{print $1,$2,$3,$4}}' | "
			"bedToIgv -name -slop 200 -path {params.path} -i stdin "
			"> {output.igv_script} 2>> {log}"

rule multi_sample_classification:
	input:
		sample_peaks = get_donor_samples
	output:
		donor_peaks = "all_{donor}_classified_peaks.bed",
		multi_peaks = "{donor}_multi_peaks.bed",
	params:
		num_samples = lambda wc, input: len(input.sample_peaks),
		comparison = get_comparison_string
	resources:
		runtime = get_min_runtime,
		mem_mb = get_filter_and_classify_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/multi_sample_classification/{donor}.log"
	group: "classify"
	shell:
		"python {workflow.basedir}/scripts/multi_classify.py "
			"-d {wildcards.donor} "
			"-f {input.sample_peaks} "
			"-p {output.donor_peaks} "
			"-m {output.multi_peaks} "
			"-c {params.comparison} "
			"2> {log}"
