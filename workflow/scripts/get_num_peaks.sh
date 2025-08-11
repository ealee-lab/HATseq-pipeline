#!/bin/bash

# Calculate the number of peaks filtered out at each step.

peak_file=$1 # Path to input peaks intersected with ground truth set (GTS)
noerr_file=$2 # Path to input peaks intersected with GTS and error-prone filtered

# echo -e "\nBEFORE ERROR-PRONE FILTER"
# echo "Intersecting peaks: "
# awk '{if ($15 > 0) print $4,$7,$8}' "${peak_file}" | sort | uniq | awk '{print $2,$3}' | sort | uniq -c

echo -e "\nAFTER ERROR-PRONE FILTER"
echo "All peaks: "
awk '{print $4,$7,$8}' "${noerr_file}" | sort | uniq | awk '{print $2,$3}' | sort | uniq -c

echo "Intersecting peaks: "
awk '{if ($15 > 0) print $4,$7,$8}' "${noerr_file}" | sort | uniq | awk '{print $2,$3}' | sort | uniq -c
