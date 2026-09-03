#!/usr/bin/env bash
 
exec 2> "${snakemake_log[0]}"

echo -e '#chrm\tstart\tend\tpeak_name\tnum_subpeaks;max_depth;length_subpeaks;depth_subpeaks\tstrand' \
	> "${snakemake_output[peaks]}" 

# if [[ "${snakemake_params[library]}" == "bulk" ]]; then dist=10; else dist=0; fi

cat <( genomeCoverageBed -ibam "${snakemake_input[peak_sorted_bam]}" -bg -ignoreD -strand + | \
	awk '{OFS="\t"; print $1,$2,$3,$4,".","+",$3-$2}' ) \
	<( genomeCoverageBed -ibam "${snakemake_input[peak_sorted_bam]}" -bg -ignoreD -strand - | \
	awk '{OFS="\t"; print $1,$2,$3,$4,".","-",$3-$2}' ) | \
	sort -k1,1 -k2,2n | \
	mergeBed -d 1 -i stdin -s -c 1,4,6,7,4 -o count,collapse,distinct,collapse,max | 
	grep -v -E "chrEBV|chrM|chrUn|random" | \
	sort -k1,1V -k2,2n | \
	awk -v sample="${snakemake_wildcards[sample]}" \
		'{ OFS="\t"; \
		if($6 == "+") \
			{print $1,$2,$3,sample"-plus-peak-"NR,$4";"$8";"$7";"$5,$6} \
		else \
			{print $1,$2,$3,sample"-minus-peak-"NR,$4";"$8";"$7";"$5,$6} }' \
	>> "${snakemake_output[peaks]}" 

# touch "${snakemake_output[peak_readID_list]}"

# while IFS= read -r line; do
# 	region=$(awk -F '\t' '{print $1":"$2"-"$3}' <<< "$line")
# 	peak=$(awk '{print $4}' <<< "$line")
# 	option=$(awk '{if($6=="+") {print "-F 16"} else {print "-f 16"}}' <<< "$line")

# 	samtools view "$option" "${snakemake_input[peak_sorted_bam]}" "$region" | \
# 		awk -v peak="$peak" '{OFS="\t"; print peak,$1}' >> "${snakemake_output[peak_readID_list]}"
# done < <(tail -n +2 "${snakemake_output[peaks]}")

bedtools getfasta -fi "${snakemake_input[ref_genome]}" -bed "${snakemake_output[peaks]}" -name \
	> "${snakemake_output[peak_seq]}" 

# awk '{ OFS="\t"; \
# 	split($5,a,";"); split(a[4],b,","); split(a[3],c,","); \
# 	sum=0; stop=1; \
# 	if($6 == "+") \
# 		{ for(i=length(b);i>0; i=i-1) \
# 			{ if(b[i] < a[2]) \
# 				{sum=sum+(stop*c[i])} \
# 			else \
# 				{stop=0} } } \
# 	else \
# 		{ for(i=1;i<length(b)+1; i=i+1) \
# 			{ if(b[i] < a[2]) \
# 				{sum=sum+(stop*c[i])} \
# 			else \
# 				{stop=0} } } print $4,sum }' "${snakemake_output[peaks]}" \
# 	> "${snakemake_output[max_depth]}" 
