import os
import re, glob
import argparse
import pandas as pd

def main(workdir, outfile, proportion):
	# Load necessary files
	all_peaks = pd.read_csv(
		f"{workdir}/all_putative_peaks.bed", sep="\t", usecols=[0,1,2,3,4,5,6])
	all_peaks["sample"] = all_peaks["peak"].str.extract(r'(.*)-(?:plus|minus)-peak-[0-9]+')
	samples = all_peaks["sample"].sort_values().unique()

	knr_peaks = glob.glob(f"{workdir}/*_KNR_peaks.bed")

	# Get number of KNRs per donor
	donors = []
	num_knrs = {}

	for bed in knr_peaks:
		donor = re.search(r'(.*)_KNR_peaks.bed', os.path.basename(bed)).group(1)
		donors.append(donor)

		with open(bed, "r") as f:
			num_knr = sum([1 for _ in f]) - 1
		num_knrs[donor] = num_knr
	
	donor_str = "|".join(donors)

	# Get sample-specific KNR recovery rates
	recall_rates = {}

	for sample in samples:
		sample_file = glob.glob(f"{workdir}/{sample}/*classified_peaks.bed")[0]
		sample_df = pd.read_csv(sample_file, sep="\t")

		donor = re.search(rf'{donor_str}', os.path.basename(sample_file)).group(0)
		sample_knr = len(
			sample_df[(sample_df["classification"] == "KNR") & (sample_df["filter"] == "PASS")])
		
		recall_rates[sample] = sample_knr / num_knrs[donor]

	# Calculate somatic L1 retrotransposition rates per sample
	count_df = all_peaks.groupby(["sample","classification"]).agg({
		'peak': lambda x: len(list(set(x))),
		'RPM': 'sum'
	}).reset_index()
	
	ratefile = open(outfile, 'w')
	
	for sample in samples:
		knr_copy = count_df.loc[
			(count_df["sample"] == sample) & (count_df["classification"] == "KNR"), "peak"].iloc[0]
		knr_read = count_df.loc[
			(count_df["sample"] == sample) & (count_df["classification"] == "KNR"), "RPM"].iloc[0]
		
		private_read = count_df.loc[
			(count_df["sample"] == sample) & (count_df["classification"] == "SOM_private"), "RPM"].iloc[0]
		clonal_read = count_df.loc[
			(count_df["sample"] == sample) & (count_df["classification"] == "SOM_clonal"), "RPM"].iloc[0]
		som_read = private_read + clonal_read
		
		# Adjust by unmasked proportion of the genome and sample-specific recovery rate
		som_copy = ((knr_copy * (som_read / knr_read)) / proportion) / recall_rates[sample]
		ratefile.write(sample + "\t" + str(som_copy) + "\n")
	
	ratefile.close()
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
        description="A script to calculate the somatic L1 retrotransposition rate per sample.")
	parser.add_argument("-d", "--directory", 
                        help="Path to Snakemake working directory with results",
                        required=True)
	parser.add_argument("-c", "--constant", 
                        help="Proportion of the genome that is not masked",
                        type=float, required=True)
	parser.add_argument("-o", "--outfile", 
                        help="Path to output file",
                        required=True)
	args = parser.parse_args()

	main(args.directory, args.outfile, args.constant)
