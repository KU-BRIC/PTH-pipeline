####################################################################################################################
####################################################################################################################
# Perform QC for sequence data with FASTQuick - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 29 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=70GB
#PBS -l walltime=200:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
Funcotator_DB="/home/projects/cu_10184/people/haikon/Reference/Funcotator/funcotator_dataSources.v1.7.20200521s"

# FASTQuick:
fastquick_reference_dir="/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37"
fastquick_resource_dir="/home/projects/cu_10184/people/haikon/Software/FASTQuick/resource"
target_fastquick=/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37/${target_name}.bed

####################################################################################################################
# Software tools:
FASTQuick="/home/projects/cu_10184/people/haikon/Software/FASTQuick/bin/FASTQuick.sh"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# FASTQuick:
# Change to working directory.
cd ${temp_dir}

# Create a folder to save index files for the sample.
mkdir ${batch_dir}/QC/FASTQuick/Lock/${sample}

# Create a folder to save the results for the sample.
mkdir -p ${batch_dir}/QC/FASTQuick/Result/${sample}/tmp

# Run FASTQuick.
${FASTQuick}  \
  --steps All  \
  --fastq_1 ${fq_dir}/${sample}_R1.fq.gz  \
  --fastq_2 ${fq_dir}/${sample}_R2.fq.gz  \
  --output ${batch_dir}/QC/FASTQuick/Result/${sample}/${sample}  \
  --workingDir ${batch_dir}/QC/FASTQuick/Result/${sample}/tmp  \
  --reference ${fastquick_reference_dir}/hs37d5.fa  \
  --targetRegion ${target_fastquick}  \
  --dbSNP ${fastquick_reference_dir}/dbsnp132_20101103.vcf.gz  \
  --callableRegion ${fastquick_reference_dir}/20141020.strict_mask.whole_genome.bed  \
  --index ${batch_dir}/QC/FASTQuick/Lock/${sample}/index  \
  --candidateVCF ${fastquick_resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --SVDPrefix ${fastquick_resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --nThread ${n_thread}

####################################################################################################################
####################################################################################################################
