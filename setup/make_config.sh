#!/bin/bash 


#Suggested command line use of this script:



sample_name=$1
fastq1=$2
fastq2=$3
library=$4
index=$5 



output_dir=/lab-share/Gene-Lee-ANR-e2/shayna/data/output//lab-share/Gene-Lee-ANR-e2/shayna/data/output/${sample_name}
tmp_dir=/temp_work/ch244189/HATseq-pipeline/tmp/${sample_name}
te_name="L1HS"
cores=8
fastq_screen_config=/lab-share/Gene-Lee-e2/Public/home/shayna/fastq_screen.conf
bwa_ref_genome=/lab-share/Gene-Lee-e2/Public/home/shayna/ref_hg38/bwa/Homo_sapiens_assembly38.fasta
STAR_genome_dir=/lab-share/Gene-Lee-e2/Public/home/shayna/ref_hg38/STAR_35_genomeDir
reference_dir=/lab-share/Gene-Lee-e2/Public/home/shayna/HATseq-pipeline/reference
min_to_pass=1
num_cells=1

echo "output_dir: \"${output_dir}\"" > ../config/config.${library}.${sample_name}
echo "tmp_dir: \"${tmp_dir}\"" >> ../config/config.${library}.${sample_name}
echo "sample_name: \"${sample_name}\"" >> ../config/config.${library}.${sample_name}
echo "te_name: \"${te_name}\"" >> ../config/config.${library}.${sample_name}
echo "index: \"${index}\"" >> ../config/config.${library}.${sample_name}
echo "cores: \"${cores}\"" >> ../config/config.${library}.${sample_name}
echo "fastq_screen_config: \"${fastq_screen_config}\"" >> ../config/config.${library}.${sample_name}
echo "bwa_ref_genome: \"${bwa_ref_genome}\"" >> ../config/config.${library}.${sample_name}
echo "STAR_genome_dir: \"${STAR_genome_dir}\"" >> ../config/config.${library}.${sample_name}
echo "reference_dir: \"${reference_dir}\"" >> ../config/config.${library}.${sample_name}
echo "min_to_pass: \"${min_to_pass}\"" >> ../config/config.${library}.${sample_name}
echo "num_cells: \"${num_cells}\"" >> ../config/config.${library}.${sample_name}
echo "library: \"${library}\"" >> ../config/config.${library}.${sample_name}
