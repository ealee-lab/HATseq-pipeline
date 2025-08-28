import pandas as pd 

def get_min_runtime(wildcards):
	return "10s"

def get_med_runtime(wildcards):
	return 6 

def get_min_mem_mb(wildcards):
	return 100

def get_qc_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 500) + 2, 30)
	except WorkflowError: # occurs during dry run
		minutes = 30
	return minutes

def get_preprocessing_runtime(wildcards, input):
	minutes = min(int(input.size_mb / 400) + 2, 60)
	return minutes

def get_alignment_runtime(wildcards, input):
	minutes = min(int((input.size_mb - 3000) / 33) + 6, 90)
	return minutes