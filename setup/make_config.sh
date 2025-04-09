#!/bin/bash 


#Suggested command line use of this script:



sample_name=$1
fastq1=$2
fastq2=$3
library=$4 #bulk, micro, single
index=$5 


seq_plat="NovaseqS4"
seq_cent="Psomagen" 
seq_plat=$6
seq_cent=$7



output_dir=/path/to/desired/output/location/${sample_name}
tmp_dir=/path/to/temp/space/${sample_name}
te_name="L1HS"
cores=8 
bwa_ref_genome=/path/to/reference/genome/fasta/Homo_sapiens_assembly38.fasta #Ensure the index and annotation files are present in the same directory  
reference_dir=/path/to/HATseq-pipeline/reference
min_to_pass=1
num_cells=1

echo "output_dir: \"${output_dir}\"" > ../config/config.${library}.${sample_name}
echo "tmp_dir: \"${tmp_dir}\"" >> ../config/config.${library}.${sample_name}
echo "sample_name: \"${sample_name}\"" >> ../config/config.${library}.${sample_name}
echo "te_name: \"${te_name}\"" >> ../config/config.${library}.${sample_name}
echo "index: \"${index}\"" >> ../config/config.${library}.${sample_name}
echo "cores: \"${cores}\"" >> ../config/config.${library}.${sample_name}
echo "raw_fastq1: \"${fastq1}\"" >> ../config/config.${library}.${sample_name}
echo "raw_fastq2: \"${fastq2}\"" >> ../config/config.${library}.${sample_name}
echo "fastq_screen_config: \"${fastq_screen_config}\"" >> ../config/config.${library}.${sample_name}
echo "bwa_ref_genome: \"${bwa_ref_genome}\"" >> ../config/config.${library}.${sample_name}
echo "STAR_genome_dir: \"${STAR_genome_dir}\"" >> ../config/config.${library}.${sample_name}
echo "reference_dir: \"${reference_dir}\"" >> ../config/config.${library}.${sample_name}
echo "min_to_pass: \"${min_to_pass}\"" >> ../config/config.${library}.${sample_name}
echo "num_cells: \"${num_cells}\"" >> ../config/config.${library}.${sample_name}
echo "library: \"${library}\"" >> ../config/config.${library}.${sample_name}
echo "sequencing_platform: \"${seq_plat}\"" >> ../config/config.${library}.${sample_name}
echo "sequencing_center: \"${seq_cent}\"" >> ../config/config.${library}.${sample_name}
