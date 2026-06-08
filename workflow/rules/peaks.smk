ref_dir = config["reference_dir"]


rule peak_calling:
	input:
		peak_sorted_bam="{sample}/{sample}_peak_sorted_bwa.bam",
		ref_genome=config["bwa_ref_genome"],
		hg38=f"{ref_dir}/human/hg38.genome",
	output:
		peaks="{sample}/{sample}_peaks.bed",
		max_depth="{sample}/{sample}_max_depth_distance_to_boundary.txt",
		# breakpoint_fa=temp("{sample}/{sample}_breakpoint.fa"),
		peak_seq="{sample}/{sample}_peak_sequence.fa",
	params:
		library=config["library"],
	resources:
		runtime=get_peak_calling_runtime,
		mem_mb=get_peak_calling_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/peak_calling/{sample}.log",
	group:
		"peaks"
	script:
		"../scripts/peak_calling.sh"


rule count_templates:
	input:
		peak_sorted_bam="{sample}/{sample}_peak_sorted_bwa.bam",
		peaks="{sample}/{sample}_peaks.bed",
	output:
		unique_start_positions="{sample}/{sample}_unique_start_positions.bed",
		peak_readID_list="{sample}/{sample}_peak_readID_list.txt",
	resources:
		runtime=get_count_templates_runtime,
		mem_mb=get_count_templates_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/count_templates/{sample}.log",
	group:
		"peaks"
	script:
		"../scripts/count_templates.py"


rule nearby_peak:
	input:
		peaks="{sample}/{sample}_peaks.bed",
	output:
		nearby_peaks="{sample}/{sample}_nearby_peaks.tsv",
	resources:
		runtime=get_min_runtime,
		mem_mb=get_min_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/nearby_peak/{sample}.log",
	group:
		"peaks"
	shell:
		'mergeBed -d 1000 -i {input.peaks} -c 4,5 -o collapse,collapse -delim "=" | '
			"awk '$4 ~ \"=\"' | "
			'awk \'{{ OFS="\\t"; split($5,a,"="); split($4,b,"="); max=0; stop=0; '
			"for(i=0; i<length(a); i++) "
			'{{split(a[i],d,";"); if(d[2] > max) {{max = d[2]; stop=i}} }} '
			"print b[stop],$4 }}' "
			"> {output.nearby_peaks} 2> {log} \n"


rule intersect_regions:
	input:
		peaks="{sample}/{sample}_peaks.bed",
		satellites=f"{ref_dir}/RepeatMasker/hg38.repeatmasker.Satellite.bed",
		segdups=f"{ref_dir}/SegDup/hg38.genomicSuperDups.v37.chr.bed",
		homopolymers=f"{ref_dir}/human/hg38.hg19.homo8.chr.bed",
	output:
		satellite_intersect="{sample}/{sample}_satellite_intersect.bed",
		segdup_intersect="{sample}/{sample}_segdup_intersect.bed",
		homopolymer_intersect="{sample}/{sample}_homopolymer_intersect.bed",
	resources:
		runtime=get_min_runtime,
		mem_mb=get_min_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/intersect_regions/{sample}.log",
	group:
		"peaks"
	script:
		"../scripts/intersect_regions.sh"


rule intersect_databases:
	input:
		peaks="{sample}/{sample}_peaks.bed",
		hg38=f"{ref_dir}/human/hg38.genome",
		repeat_masker=f"{ref_dir}/RepeatMasker/hg38.repeatmasker.L1.bed",
		Evrony_KR=f"{ref_dir}/RepeatMasker/hg38.Evrony_KR_960.liftover.bed",
		# satellites=f"{ref_dir}/RepeatMasker/hg38.repeatmasker.Satellite.bed",
		i1kgp=f"{ref_dir}/1kgp/ALL_MELT_ME_1000G_HC_20190901.AF.bed",
		gnomad=f"{ref_dir}/gnomAD-SV/gnomad.v4.1.ME.sites.bed",
		nyuwa=f"{ref_dir}/nyuwa/MEI.GRCh38.HMEIDv1.1.final.bed",
		xtea=f"{ref_dir}/xTea/xTea_Borges-Monroy2021.accessioned.bed",
		hgsvc3=f"{ref_dir}/HGSVC3/MEI_Callset_GRCh38.ALL.20241211.bed",
		melt_lra=f"{ref_dir}/HGSVC3/Ortho_MEI_GRCh38.ALL.20241211.bed",
		ont=f"{ref_dir}/1019_ONT/1019_ONT_Schloissnig2025.bed",
	output:
		intersect_annotated="{sample}/{sample}_intersect_annotated.bed",
	resources:
		runtime=get_min_runtime,
		mem_mb=get_intersect_databases_mem_mb,
	conda:
		"../envs/HATseq.yml"
	log:
		"logs/intersect_databases/{sample}.log",
	group:
		"peaks"
	script:
		"../scripts/intersect_databases.sh"
