import pandas as pd

#### Manipulate metadata ####
def set_sample_names(samples, library):
	samples.insert(0, "sample_name", None) 

	if library == "bulk":
		samples.loc[samples["tissue"] == "-NA-", "sample_name"] = samples["donor"]
		samples.loc[samples["tissue"] != "-NA-", "sample_name"] = samples["donor"] + "-" + samples["tissue"]
	else:
		samples.loc[samples["cell"] == "-NA-", "sample_name"] = samples["donor"]
		samples.loc[samples["cell"] != "-NA-", "sample_name"] = samples["donor"] + "-" + samples["cell"]
	return

def get_donor_samples(wildcards, samples):
	names = samples[samples["donor"] == wildcards.donor]["sample_name"]
	return names

#### Define memory limits for each rule ####
def get_min_mem_mb(wildcards, attempt):
	mb = 120 + (100 * (attempt - 1))
	return mb

def get_qc_mem_mb(wildcards, attempt):
	mb = 2000 + (1000 * (attempt - 1))
	return mb

def get_preprocessing_mem_mb(wildcards, input, attempt):
	mb = max(0.75 * input.size_mb + (4000 * (attempt - 1)), 2000)
	return mb

def get_alignment_mem_mb(wildcards, attempt):
	mb = 12000 + (2000 * (attempt - 1))
	return mb

def get_Ntag_mem_mb(wildcards, input, attempt):
	mb = max((2.5 * input.size_mb - 5000) + (2500 * (attempt - 1)), 2000)
	return mb

def get_gmotif_mem_mb(wildcards, input, attempt):
	mb = max((2.75 * input.size_mb - 2000) + (2000 * (attempt - 1)), 1250)
	return mb

def get_extract_clipping_mem_mb(wildcards, attempt):
	mb = 150 + (100 * (attempt - 1))
	return mb

def get_peak_calling_mem_mb(wildcards, attempt):
	mb = 4000 + (2000 * (attempt - 1))
	return mb

def get_count_templates_mem_mb(wildcards, input, attempt):
	mb = max((5 * input.size_mb - 600) + (1000 * (attempt - 1)), 1000)
	return mb

def get_intersect_databases_mem_mb(wildcards, attempt):
	# TODO: Memory spikes in this rule, resolve if possible
	mb = 2000 + (1000 * (attempt - 1))
	return mb

def get_create_big_table_mem_mb(wildcards, input, attempt):
	mb = max((2.75 * input.size_mb) + (2000 * (attempt - 1)), 2000)
	return mb

def get_filter_and_classify_mem_mb(wildcards, attempt):
	mb = 600 + (300 * (attempt - 1))
	return mb

def get_map_transduction_fasta_mem_mb(wildcards, attempt):
	mb = 2000 + (2000 * (attempt - 1))
	return mb

def get_transduction_big_table_mem_mb(wildcards, input, attempt):
	mb = 8000 + (1000 * (attempt - 1))
	return mb

def get_call_transductions_mem_mb(wildcards, attempt):
	mb = 600 + (150 * (attempt - 1))
	return mb

#### Define runtime limits for each rule ####
def get_min_runtime(wildcards, attempt):
	minutes = 1 + (1 * (attempt - 1))
	return minutes
	
def get_qc_runtime(wildcards, input, attempt):
	minutes = int(input.size_mb / 800) + 10 + (10 * (attempt - 1))
	return minutes

def get_preprocessing_runtime(wildcards, input, attempt):
	try:
		minutes = int(input.size_mb / 400) + 12 + (15 * (attempt - 1))
	except FileNotFoundError: # occurs during dry run
		minutes = 60
	return minutes

def get_alignment_runtime(wildcards, input, attempt):
	try:
		minutes = int(input.size_mb / 80) + 10 + (20 * (attempt - 1))
	except FileNotFoundError:
		minutes = 90
	return minutes

def get_Ntag_runtime(wildcards, input, attempt):
	try:
		minutes = min(int(input.size_mb / 333) - 4.5 + (5 * (attempt - 1)), 20)
	except FileNotFoundError:
		minutes = 20
	return minutes

def get_gmotif_runtime(wildcards, input, attempt):
	try:
		minutes = int(input.size_mb / 333) + 2 + (5 * (attempt - 1))
	except FileNotFoundError:
		minutes = 15
	return minutes

def get_junction_spanning_runtime(wildcards, input, attempt):
	try:
		minutes = int(input.size_mb / 500) + 2 + (2 * (attempt - 1))
	except FileNotFoundError:
		minutes = 10
	return minutes

def get_extract_clipping_runtime(wildcards, input, attempt):
	try:
		minutes = int(input.size_mb / 400) + 2 + (2 * (attempt - 1))
	except FileNotFoundError:
		minutes = 10
	return minutes

def get_peak_calling_runtime(wildcards, input, attempt):
	try:
		minutes = int(input.size_mb / 300) + 2 + (2 * (attempt - 1))
	except FileNotFoundError:
		minutes = 10
	return minutes

def get_count_templates_runtime(wildcards, input, attempt):
	try:
		minutes = max(int(input.size_mb / 375) + 1 + (2 * (attempt - 1)), 1)
	except FileNotFoundError:
		minutes = 10
	return minutes

def get_create_big_table_runtime(wildcards, input, attempt):
	try:
		minutes = max(int(input.size_mb / 700) + 1 + (3 * (attempt - 1)), 2)
	except FileNotFoundError:
		minutes = 15
	return minutes

def get_map_transduction_fasta_runtime(wildcards, input, attempt):
	minutes = 2 + (1 * (attempt - 1))
	return minutes

def get_transduction_big_table_runtime(wildcards, input, attempt):
	minutes = 4 + (1 * (attempt - 1))
	return minutes
