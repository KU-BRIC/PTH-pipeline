####################################################################################################################
####################################################################################################################
# Perform sequence and sample quanlity control - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 6 May 2021
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
reference_dir="/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37"
resource_dir="/home/projects/cu_10184/people/haikon/Software/FASTQuick/resource"
index_dir="/home/projects/cu_10184/people/haikon/Software/FASTQuick/index"

# Software tools:
verifyBamID="/home/projects/cu_10184/people/haikon/Software/VerifyBamID/verifyBamID"

# VCF and BAM file
vcf=/home/projects/cu_10184/projects/PTH/QC/Lock/FingerPrinting/VerifyBamID/Genotype/${batch2}/${sample2}.vcf
bam=${bam_dir}/${sample}.bam

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Create a folder to save the results.
rm -rf ${QC_dir}/${sample}_vs_${sample2}
mkdir -p ${QC_dir}/${sample}_vs_${sample2}/tmp

rm -rf ${index_dir}/${sample}_vs_${sample2}
mkdir -p ${index_dir}/${sample}_vs_${sample2}

${verifyBamID} --vcf ${vcf} --bam ${bam} --out ${QC_dir}/${sample}_vs_${sample2}/${sample}

####################################################################################################################
####################################################################################################################

