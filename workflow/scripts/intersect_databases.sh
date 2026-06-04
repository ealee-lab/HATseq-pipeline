#!/usr/bin/env bash
 
exec 2> "${snakemake_log[0]}"

echo -e '#chr\tstart\tend\t'\
'peak_name\tnum_subpeaks;max_depth;length_subpeaks;depth_subpeaks\tstrand\t'\
'RepeatMasker\tEvrony_KR\t'\
'1000_Genomes_Project\tgnomAD\tNyuWa\txTea\tHGSVC3\tHGSVC3-MELT-LRA\t1019_ONT' \
> "${snakemake_output[intersect_annotated]}" 

intersectBed -S -wao -a "${snakemake_input[peaks]}" \
	-b <(slopBed -g "${snakemake_input[hg38]}" -l 0 -r 500 -s -i "${snakemake_input[repeat_masker]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$10";"$11";"$12}' | \
	groupBy -grp 1,2,3,4,5,6 -c 7 -o collapse | \
	intersectBed -S -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -l 0 -r 500 -s -i "${snakemake_input[Evrony_KR]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$11";"$12";"$13}' | \
	groupBy -grp 1,2,3,4,5,6,7 -c 8 -o collapse | \
	intersectBed -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -b 500 -i "${snakemake_input[i1kgp]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$12";"$13";"$14}' | \
	groupBy -grp 1,2,3,4,5,6,7,8 -c 9 -o collapse | \
	intersectBed -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -b 500 -i "${snakemake_input[gnomad]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13";"$14}' | \
	groupBy -grp 1,2,3,4,5,6,7,8,9 -c 10 -o collapse | \
	intersectBed -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -b 500 -i "${snakemake_input[nyuwa]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14";"$15}' | \
	groupBy -grp 1,2,3,4,5,6,7,8,9,10 -c 11 -o collapse | \
	intersectBed -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -b 500 -i "${snakemake_input[xtea]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$15";"$16";"$17}' | \
	groupBy -grp 1,2,3,4,5,6,7,8,9,10,11 -c 12 -o collapse | \
	intersectBed -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -b 500 -i "${snakemake_input[hgsvc3]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16";"$17}' | \
	groupBy -grp 1,2,3,4,5,6,7,8,9,10,11,12 -c 13 -o collapse | \
	intersectBed -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -b 500 -i "${snakemake_input[melt_lra]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$17";"$18}' | \
	groupBy -grp 1,2,3,4,5,6,7,8,9,10,11,12,13 -c 14 -o collapse | \
	intersectBed -wao -a stdin \
	-b <(slopBed -g "${snakemake_input[hg38]}" -b 500 -i "${snakemake_input[ont]}") | \
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$18}' | \
	groupBy -grp 1,2,3,4,5,6,7,8,9,10,11,12,13,14 -c 15 -o collapse \
	>> "${snakemake_output[intersect_annotated]}"
