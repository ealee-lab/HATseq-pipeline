import pandas as pd 

def get_short_runtime(wildcards):
	try:
		minutes = min(int(input.size_mb / 75) + 5, 60)
	except:
		minutes = 60
	return minutes

def get_medium_runtime(wildcards):
	try:
		minutes = min(int(input.size_mb / 100) + 5, 90)
	except:
		minutes = 90
	return minutes

def get_long_runtime(wildcards):
	try:
		minutes = min(int(input.size_mb / 50) + 10, 150)
	except:
		minutes = 150
	return minutes