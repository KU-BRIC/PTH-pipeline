####################################################################################################################
####################################################################################################################
# Perform sample quanlity control with FASTQuick - job.
# Author: Balthasar Schlotmann
# Last Modified: 10 June 2021
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

# Software tools:
FASTQuick="/home/projects/cu_10184/people/haikon/Software/FASTQuick/bin/FASTQuick.sh"

export PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin/:$PATH

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Create a folder to save index files.
mkdir -p ${lock_dir}/${sample}

# Create a folder to save the results.
mkdir -p ${res_dir}/${sample}/tmp

####################################################################################################################
# Run FASTQuick.
${FASTQuick}  \
  --steps All  \
  --fastq_1 ${fq_dir}/${sample}_R1.fq.gz  \
  --fastq_2 ${fq_dir}/${sample}_R2.fq.gz  \
  --output ${res_dir}/${sample}/${sample}  \
  --workingDir ${res_dir}/${sample}/tmp  \
  --reference ${reference_dir}/hs37d5.fa  \
  --targetRegion ${target}  \
  --dbSNP ${reference_dir}/dbsnp132_20101103.vcf.gz  \
  --callableRegion ${reference_dir}/20141020.strict_mask.whole_genome.bed  \
  --index ${lock_dir}/${sample}/index  \
  --candidateVCF ${resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --SVDPrefix ${resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --nThread 8

####################################################################################################################
# Convert html to pdf.
# pandoc ${res_dir}/${sample}/${sample}.FinalReport.html -t latex -o ${res_dir}/${sample}/${sample}.FinalReport.pdf

####################################################################################################################
####################################################################################################################

