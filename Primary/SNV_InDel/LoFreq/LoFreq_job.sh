####################################################################################################################
####################################################################################################################
# Run LoFreq and annotate - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 14 July 2021
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

# Panel target:
target_chr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${batch_dir}/temp

####################################################################################################################
# Call SNV_InDel variants and annotate.
####################################################################################################################
# LoFreq:
################################################
# Set parameters.
min_bq=20
min_alt_bq=20
min_mq=20
sig=0.01
min_cov=1

# Run LoFreq.
lofreq call-parallel --pp-threads ${n_thread} --call-indels -f $hg -l ${target_chr}  \
  -q ${min_bq} -Q ${min_alt_bq} -m ${min_mq} -a ${sig} -C ${min_cov}  \
  -o ${batch_dir}/Lock/SNV_InDel/LoFreq/vcf/${sample}.vcf ${batch_dir}/Lock/BAM/${sample}.bam

# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${batch_dir}/Lock/SNV_InDel/LoFreq/vcf/${sample}.vcf -O ${batch_dir}/Lock/SNV_InDel/LoFreq/maf/${sample}.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Fill_DP_AF_FuncotatorMAF_nonSNVer_job.R ${batch_dir}/Lock/SNV_InDel/LoFreq ${sample} 'LoFreq'

####################################################################################################################
####################################################################################################################
