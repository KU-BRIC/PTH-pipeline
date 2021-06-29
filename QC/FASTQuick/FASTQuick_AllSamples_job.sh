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
# all samples

sample="all_samples"

####################################################################################################################
####################################################################################################################
# Reference databases:
reference_dir="/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37"
resource_dir="/home/projects/cu_10184/people/haikon/Software/FASTQuick/resource"
index_dir="/home/projects/cu_10184/people/haikon/Software/FASTQuick/index"

# Software tools:
FASTQuick="/home/projects/cu_10184/people/haikon/Software/FASTQuick/bin/FASTQuick.sh"
# FASTQuick="/home/projects/cu_10145/people/balsch/PTH/FASTQuick/FASTQuick/bin/FASTQuick.sh"

export PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin/:$PATH

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Create a folder to save the results.
rm -rf ${QC_dir}/${sample}
mkdir -p ${QC_dir}/${sample}/tmp

rm -rf ${index_dir}/${sample}
mkdir -p ${index_dir}/${sample}

${FASTQuick}  \
  --steps All  \
  --fastqList ${QC_dir}/fastqfile.tsv \
  --output ${QC_dir}/${sample}/${sample}  \
  --workingDir ${QC_dir}/${sample}/tmp  \
  --reference ${reference_dir}/hs37d5.fa  \
  --targetRegion ${target}  \
  --dbSNP ${reference_dir}/dbsnp132_20101103.vcf.gz  \
  --callableRegion ${reference_dir}/20141020.strict_mask.whole_genome.bed  \
  --index ${index_dir}/${sample}/index  \
  --candidateVCF ${resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --SVDPrefix ${resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --nThread 8

####################################################################################################################
# convert to pdf

#pandoc ${QC_dir}/${sample}/${sample}.FinalReport.html -t latex -o ${QC_dir}/${sample}/${sample}.FinalReport.pdf

####################################################################################################################
####################################################################################################################

