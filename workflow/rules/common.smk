import pandas as pd 

def get_min_mem_mb(wildcards):
	return 100

def get_qc_mem_mb(wildcards):
	return 2000

def get_preprocessing_mem_mb(wildcards):
	return 8000

def get_alignment_mem_mb(wildcards):
	return 12000

def get_Ntag_mem_mb(wildcards, input):
	# TODO: Update mem
	# mb = max(1.75 * input.size_mb, 2000)
	return 2400

def get_gmotif_mem_mb(wildcards, input):
	mb = max(2.5 * input.size_mb, 1000)
	return mb

# def get_junction_spanning_mem_mb(wildcards):
# 	return 100

def get_extract_clipping_mem_mb(wildcards, input):
	return 150

def get_peak_calling_mem_mb(wildcards, input):
	mb = max(2 * input.size_mb, 2000)
	return mb

def get_count_templates_mem_mb(wildcards, input):
	mb = max(6 * input.size_mb, 2000)
	return mb

def get_intersect_databases_mem_mb(wildcards, input):
	# TODO: Memory spikes in this rule, resolve if possible
	return 4000

def get_create_big_table_mem_mb(wildcards, input):
	mb = max(4 * input.size_mb, 4000)
	return mb

def get_filter_and_classify_mem_mb(wildcards):
	return 1000

def get_map_transduction_fasta_mem_mb(wildcards):
	return 2000

def get_transduction_big_table_mem_mb(wildcards, input):
	mb = max(4 * input.size_mb, 2000)
	return mb

def get_call_transductions_mem_mb(wildcards):
	return 600

def get_min_runtime(wildcards):
	return 1
	
def get_qc_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 500) + 5, 30)
	except FileNotFoundError: # occurs during dry run
		minutes = 30
	return minutes

def get_preprocessing_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 400) + 5, 60)
	except FileNotFoundError:
		minutes = 60
	return minutes

def get_alignment_runtime(wildcards, input):
	try:
		minutes = min(int((input.size_mb - 3000) / 33) + 5, 90)
	except FileNotFoundError:
		minutes = 90
	return minutes

def get_Ntag_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 1000) + 2, 20)
	except FileNotFoundError:
		minutes = 20
	return minutes

def get_gmotif_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 333) + 1, 12)
	except FileNotFoundError:
		minutes = 12
	return minutes

def get_junction_spanning_runtime(wildcards):
	try:
		minutes = min(int(input.size_mb / 500) + 2, 6)
	except FileNotFoundError:
		minutes = 6
	return minutes

def get_extract_clipping_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 400) + 2, 10)
	except FileNotFoundError:
		minutes = 10
	return minutes

def get_peak_calling_runtime(wildcards):
	return 6

def get_count_templates_runtime(wildcards):
	return 6

def get_create_big_table_runtime(wildcards):
	# TODO: Make dynamic
	return 15

def get_filter_and_classify_runtime(wildcards):
	return 2

def get_map_transduction_fasta_runtime(wildcards):
	# TODO: Make dynamic
	return 60

def get_transduction_big_table_runtime(wildcards):
	# TODO: Make dynamic
	return 15
