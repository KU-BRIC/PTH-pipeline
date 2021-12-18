####################################################################################################################
####################################################################################################################
# Fingerprinting with Somalier - job.
# Author: Haiying Kong
# Last Modified: 29 October 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=10GB
#PBS -l walltime=20:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
sites="/home/projects/cu_10184/people/haikon/Software/somalier-0.2.13/sites.hg38.vcf.gz"
target="/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/Focused_myeloid_panel-All_target_segments_covered_by_probes-TE-93310852_hg38_v2_190722165759.bed"

# Software tools:
Somalier="/home/projects/cu_10184/people/haikon/Software/somalier-0.2.13/somalier"

# Filter bam file for small panel.
bedtools intersect -abam ${bam} -b ${target} > ${temp_dir}/${sam}.bam
samtools view -H ${bam} | grep "^@RG" > ${temp_dir}/${sam}.RG
IFS=$'\t' read -d '' -r -a rg < ${temp_dir}/${sam}.RG
rg=$(IFS=' ' ; echo "${rg[*]}")
rg="${rg// /\t}"
samtools addreplacerg -r $rg -o ${temp_dir}/${sam}_1.bam ${temp_dir}/${sam}.bam
mv ${temp_dir}/${sam}_1.bam ${temp_dir}/${sam}.bam
samtools index -b ${temp_dir}/${sam}.bam

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Filter bam file for small panel.
bedtools intersect -abam ${bam} -b ${target} > ${temp_dir}/${sam}.bam
samtools view -H ${bam} | grep "^@RG" > ${temp_dir}/${sam}.RG
IFS=$'\t' read -d '' -r -a rg < ${temp_dir}/${sam}.RG 
rg=$(IFS=' ' ; echo "${rg[*]}")
rg="${rg// /\t}"
samtools addreplacerg -r $rg -o ${temp_dir}/${sam}_1.bam ${temp_dir}/${sam}.bam
mv ${temp_dir}/${sam}_1.bam ${temp_dir}/${sam}.bam
samtools index -b ${temp_dir}/${sam}.bam

####################################################################################################################
# Somalier.
${Somalier} extract -d ${dir_name}/Lock/Somalier/Extracted --sites ${sites} -f $hg ${temp_dir}/${sam}.bam

####################################################################################################################
# Remove temp files.
rm ${temp_dir}/${sam}.*

####################################################################################################################
####################################################################################################################
