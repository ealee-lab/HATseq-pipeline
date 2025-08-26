import pandas as pd 

def get_min_runtime(wildcards):
	return "10s"

def get_med_runtime(wildcards):
	return 6 

def get_min_mem_mb(wildcards):
	return 100

def get_short_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 75) + 5, 60)
	except:
		minutes = 60
	return minutes

def get_medium_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 100) + 5, 90)
	except:
		minutes = 90
	return minutes

def get_long_runtime(wildcards, input):
	try:
		minutes = min(int(input.size_mb / 50) + 5, 90)
	except:
		minutes = 90
	return minutes