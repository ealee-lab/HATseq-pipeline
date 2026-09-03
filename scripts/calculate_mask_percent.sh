#!/bin/bash

error_prone=$1

# Run in samtools conda environment
cat <(awk -v OFS='\t' '{if($6=="+") print $0}' ../reference/all_KR.bed | \
	  grep -E 'L1Hs|L1HS|L1PA2|L1PA3|L1PA4|L1PA5' | \
	  grep -v -E 'alt|random|chrUn|chrM' | \
	  bedtools slop -s -l 0 -r 500 -i stdin -g ../reference/human/hg38.sorted.genome) \
	<(bedtools slop -b 500 -i ../reference/all_KNR.bed -g ../reference/human/hg38.sorted.genome) \
	<(awk -v OFS='\t' '{print $1,$2,$3,$4,$5,"."}' "$error_prone") | \
	sort -k1,1 -k2,2n | \
	bedtools merge -i stdin \
	> ../reference/masked_forward.bed

cat <(awk -v OFS='\t' '{if($6=="-") print $0}' ../reference/all_KR.bed | \
	  grep -E 'L1Hs|L1HS|L1PA2|L1PA3|L1PA4|L1PA5' | \
	  grep -v -E 'alt|random|chrUn|chrM' | \
	  bedtools slop -s -l 0 -r 500 -i stdin -g ../reference/human/hg38.sorted.genome) \
	<(bedtools slop -b 500 -i ../reference/all_KNR.bed -g ../reference/human/hg38.sorted.genome) \
	<(awk -v OFS='\t' '{print $1,$2,$3,$4,$5,"."}' "$error_prone") | \
	sort -k1,1 -k2,2n | \
	bedtools merge -i stdin \
	> ../reference/masked_reverse.bed

half=$(grep -v -E 'alt|random|chrUn|chrM' ../reference/human/hg38.sorted.genome | awk '{sum+=$2} END {print sum}')
total=$(($half*2))

forward=$(awk '{sum+=$3-$2} END {print sum}' ../reference/masked_forward.bed)
reverse=$(awk '{sum+=$3-$2} END {print sum}' ../reference/masked_reverse.bed)

mask=$(($forward + $reverse))
mask_pct=$(echo "scale=4; ($mask/$total)*100" | bc)

echo "Total masked percent: $mask_pct%"
