####################################################################################################################
####################################################################################################################
# Run LoFreq - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 31 May 2021
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

# Software tools:

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# LoFreq:
################################################
# Set parameters.
min_bq=6
min_alt_bq=6
min_mq=0
sig=0.01
min_cov=1

# Run LoFreq.
lofreq call-parallel --pp-threads ${n_thread} --call-indels -f $hg -l ${target_chr}  \
  -q ${min_bq} -Q ${min_alt_bq} -m ${min_mq} -a ${sig} -C ${min_cov}  \
  -o ${Lock_LoFreq_dir}/vcf/${sample}.vcf ${BAM_dir}/${sample}.bam

# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${Lock_LoFreq_dir}/vcf/${sample}.vcf -O ${Lock_LoFreq_dir}/maf/${sample}.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Fill_DP_AF_FuncotatorMAF_nonSNVer_job.R ${Lock_LoFreq_dir} ${sample} 'LoFreq'

####################################################################################################################
####################################################################################################################
