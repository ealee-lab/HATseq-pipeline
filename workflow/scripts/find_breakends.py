import sys
import pandas as pd
import numpy as np
import statistics
from natsort import natsort_key

sys.stderr = open(snakemake.log[0], "w", buffering=1)

peaks = pd.read_table(
	snakemake.input.peaks, sep="\t", skiprows=1, header=None, 
    usecols=[0,1,2,3,5], names=["chrm","start","end","peak_name","strand"])

readIDs = pd.read_table(
	snakemake.input.peak_readID_list, sep="\t", header=None, names=["peak_name","read"])
peaks = peaks.merge(readIDs, how="left", on="peak_name")
del readIDs

clips = pd.read_table(snakemake.input.softclip3p_table, sep="\t", header=None,
					  usecols=[0,1,2,3], names=["read","chrm","clippoint","strand"])
peaks = peaks.merge(clips, how="left", on=["read","chrm","strand"])
del clips

# Remove rows where read is incorrectly matched to peak
idxs = peaks[(peaks["clippoint"] >= peaks["start"]) & (peaks["clippoint"] <= peaks["end"])].index
peaks.loc[~(peaks.index.isin(idxs)), "clippoint"] = np.nan
# peaks.loc[(peaks["clippoint"] < peaks["start"]) & (peaks["strand"] == "-"), "clippoint"] = np.nan

# If peak is not 3' clipped, set clip point to 3' end of peak
peaks.loc[(peaks["clippoint"].isna()) & (peaks["strand"] == "+"), "clippoint"] = peaks["end"]
peaks.loc[(peaks["clippoint"].isna()) & (peaks["strand"] == "-"), "clippoint"] = peaks["start"]

# peaks = peaks[(peaks["clippoint"] >= peaks["start"]) & (peaks["clippoint"] <= peaks["end"])]

breakends = peaks.groupby(["chrm","peak_name","strand"]).agg({
    'clippoint': lambda x: statistics.mode(x)
}).reset_index()
del peaks

breakends["start"] = breakends["clippoint"].astype(int) - 15
breakends["end"] = breakends["clippoint"].astype(int) + 15
breakends["score"] = "."

breakends = breakends[["chrm","start","end","peak_name","score","strand"]]
breakends = breakends.sort_values(by=["chrm","start"], key=natsort_key)
breakends = breakends.rename(columns={"chrm": "#chrm"})
breakends.to_csv(snakemake.output.breakends, sep="\t", index=False, header=True)
