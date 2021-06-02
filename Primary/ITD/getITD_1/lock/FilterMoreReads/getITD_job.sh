####################################################################################################################
####################################################################################################################
# Identify ITD with getITD - job.
# Author: Haiying Kong
# Last Modified: 28 April 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=100:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
target="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/FLT3.bed"
reference="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon.txt"
anno="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv"

# Software tools:
filter_reads="/home/projects/cu_10184/projects/PTH/Code/Primary/ITD/getITD/filter_reads.py"
getITD="python3 /home/projects/cu_10184/people/haikon/Software/getITD/getitd.py"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Filter fastq manually.
mkdir ${Lock_getITD_dir}/${sample}
cd ${Lock_getITD_dir}/${sample}

# get IDs from target region plus unmapped reads, remove duplicates, because mates have same ID
{ samtools view -L ${target} ${BAM_dir}/${sample}.bam & samtools view -f 4 ${BAM_dir}/${sample}.bam; } | cut -d$'\t' -f 1 | awk '!a[$0]++' > filter_IDs.txt
gunzip -c ${fq_dir}/${sample}_R1.fq.gz | python3 ${filter_reads} filter_IDs.txt FLT3_R1.fq
gunzip -c ${fq_dir}/${sample}_R2.fq.gz | python3 ${filter_reads} filter_IDs.txt FLT3_R2.fq

####################################################################################################################
# Run getITD:
# cd ${Lock_getITD_dir}/${sample}
conda activate getITD
$getITD -reference $reference -anno $anno -infer_sense_from_alignment True -plot_coverage True -require_indel_free_primers False -min_read_length 50 -min_bqs 20 -min_read_copies 1 -filter_ins_vaf 0.001 ${sample} FLT3_R1.fq FLT3_R2.fq
conda deactivate

####################################################################################################################
# Clean the output folder.
rm *.*
mv ${sample}_getitd/* ./
rm -r ${sample}_getitd

####################################################################################################################
####################################################################################################################
