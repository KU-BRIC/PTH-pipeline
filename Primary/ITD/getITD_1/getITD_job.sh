####################################################################################################################
####################################################################################################################
# Identify ITD with modified getITD that does not exclude non-in-frame - job.
# Author: Haiying Kong
# Last Modified: 31 May 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=50:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
target_flt3="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/FLT3.bed"
getITD_reference="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon.txt"
getITD_anno="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv"

# Software tools:
getITD="python3 /home/projects/cu_10184/people/haikon/Software/getITD/getitd_1.py"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Filter bam file with target bed file and conver it to fastq.
mkdir ${Lock_getITD_dir}/$sample
cd ${Lock_getITD_dir}/$sample
bedtools intersect -abam ${BAM_dir}/${sample}.bam -b ${target_flt3} -u > FLT3.bam
# samtools view -b ${BAM_dir}/${sample}.bam chr13:28033736-28034557 > FLT3.bam
samtools sort -n FLT3.bam > tmp.bam
mv tmp.bam FLT3.bam
bedtools bamtofastq -i FLT3.bam -fq FLT3_R1.fq -fq2 FLT3_R2.fq

# Run getITD.
conda activate getITD
$getITD -reference ${getITD_reference} -anno ${getITD_anno} -infer_sense_from_alignment True -plot_coverage True -require_indel_free_primers False -min_read_length 50 -min_bqs 20 -min_read_copies 1 -filter_ins_vaf 0.001 ${sample} FLT3_R1.fq FLT3_R2.fq
conda deactivate

# Clean the output folder.
rm FLT3*
mv ${sample}_getitd/* ./
rm -r ${sample}_getitd

# Clean the output tsv.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/getITD_Clean.R ${batch} ${sample} ${Lock_getITD_dir}

####################################################################################################################
####################################################################################################################
