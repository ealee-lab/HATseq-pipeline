import re
import argparse
import pandas as pd
import pybedtools
from pybedtools import BedTool
from natsort import natsort_key

def main(filenames, outfile):
	dfs = []
	for f in filenames:
		df = pd.read_csv(f, sep="\t")
		df.columns = df.columns.str.strip("#")
		df["donor"] = re.match(r'all_(.*)_classified_peaks\.bed', f).group(1)
		dfs.append(df)
	
	peaks_df = pd.concat(dfs)
	peaks_df = peaks_df.sort_values(by=["chrm","start"])

	peaks_bt = BedTool.from_dataframe(peaks_df[peaks_df["classification"].str.contains(r'SOM|UNK')])
	merged_peaks = peaks_bt.merge(d=50, s=True, c=[4,7,9], o=["collapse","distinct","distinct"]).to_dataframe()
	merged_peaks.columns = ["chrm","start","end","peak","classification","donor"]
	merged_peaks["donor"] = merged_peaks["donor"].str.split(",").str.len()

	# Relabel UNK peaks called as SOM due to low signal in one donor
	rename_peaks = merged_peaks[(merged_peaks["classification"].str.contains("SOM")) & 
								(merged_peaks["classification"].str.contains("UNK"))]
	rename_peak_list = ",".join(rename_peaks["peak"]).split(",")
	peaks_df.loc[peaks_df["peak"].isin(rename_peak_list), "classification"] = "UNK"
	peaks_df.loc[peaks_df["peak"].isin(rename_peak_list), "filter"] = "PASS (multi-donor)"

	# Filter out low-level peaks that are recurrent across donors (likely noise)
	put_peaks = merged_peaks[merged_peaks["donor"] == 1]
	put_peak_list = ",".join(put_peaks["peak"]).split(",")

	out_df = peaks_df[
		((peaks_df["filter"].str.contains("PASS")) & 
		(~peaks_df["classification"].str.contains("SOM"))) | 
		(peaks_df["peak"].isin(put_peak_list))]

	out_df = out_df.drop(columns="donor")
	out_df = out_df.sort_values(by=["chrm","start"], key=natsort_key)
	out_df.to_csv(outfile, sep="\t", index=False, header=True)
	
	return out_df
		 
if __name__ == "__main__":
	parser = argparse.ArgumentParser(
        description="A script to filter low-level peaks recurrent across donors.")
    
	parser.add_argument("-f", "--filenames", 
                        help="List of file names of all donor-specific peaks",
                        nargs="+", 
                        required=True)
	parser.add_argument("-o", "--outfile", 
                        help="Path to output file with all putative peaks for batch",
                        required=True)
	args = parser.parse_args()

	main(args.filenames, args.outfile)
    