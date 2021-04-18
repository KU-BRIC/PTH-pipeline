####################################################################################################################
####################################################################################################################
# Identify ITD with getITD - job.
# Author: Haiying Kong
# Last Modified: 18 April 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=10:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
reference="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon.txt"
anno="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv"

# Software tools:
getITD="python3 /home/projects/cu_10184/people/haikon/Software/getITD/getitd.py"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Filter bam file with target bed file and conver it to fastq.
mkdir ${Lock_getITD_dir}/${sample}
cd ${Lock_getITD_dir}/${sample}
# bedtools intersect -abam ${BAM_dir}/${sample}.bam -b $target -u > FLT3.bam
samtools view -b ${BAM_dir}/${sample}.bam chr13:28033736-28034557 > FLT3.bam
samtools sort -n FLT3.bam > tmp.bam
mv tmp.bam FLT3.bam
bedtools bamtofastq -i FLT3.bam -fq FLT3_R1.fq -fq2 FLT3_R2.fq

####################################################################################################################
# Run getITD:
# cd ${Lock_getITD_dir}/${sample}
conda activate getITD
$getITD -reference $reference -anno $anno -infer_sense_from_alignment True -plot_coverage True -require_indel_free_primers False -min_read_length 70 -min_bqs 20 -min_read_copies 1 ${sample} FLT3_R1.fq FLT3_R2.fq
conda deactivate

####################################################################################################################
# Clean the output folder.
rm FLT3.*
mv ${sample}_getitd/* ./
rm -r ${sample}_getitd

####################################################################################################################
####################################################################################################################
