####################################################################################################################
####################################################################################################################
# Pre-process, call variants, annotate and filter variants - job.
# Author: Haiying Kong
# Last Modified: 18 July 2021
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
target_nochr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
target_chr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed
target_nopad=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Original/${target_name}.bed

####################################################################################################################
# Software tools:
picard="java -jar /home/projects/cu_10184/people/haikon/Software/picard.jar"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${batch_dir}/temp

####################################################################################################################
# Preprocessing with BWA, Picard and GATK.
####################################################################################################################
if [[ -d "/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}/Lock/BAM" ]]
then
  rm -r ${batch_dir}/Lock/BAM
  ln -s /home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}/Lock/BAM ${batch_dir}/Lock/BAM
else
  #########################################
  # Map to hg38:
  bwa mem \
    $hg \
    ${fq_dir}/${sample}_R1.fq.gz ${fq_dir}/${sample}_R2.fq.gz \
    -t ${n_thread} \
    > ${batch_dir}/Lock/BAM/${sample}_bwa.sam

  $picard FastqToSam \
    F1=${fq_dir}/${sample}_R1.fq.gz \
    F2=${fq_dir}/${sample}_R2.fq.gz \
    O=${batch_dir}/Lock/BAM/${sample}_unmapped.sam \
    PL=ILLUMINA \
    SM=${sample}

  $picard MergeBamAlignment \
    R=$hg \
    ALIGNED=${batch_dir}/Lock/BAM/${sample}_bwa.sam \
    UNMAPPED=${batch_dir}/Lock/BAM/${sample}_unmapped.sam \
    O=${batch_dir}/Lock/BAM/${sample}_mapped.bam

  rm ${batch_dir}/Lock/BAM/${sample}_bwa.sam
  rm ${batch_dir}/Lock/BAM/${sample}_unmapped.sam

  #########################################
  # Mark duplicates and sort:
  gatk MarkDuplicatesSpark \
    -I ${batch_dir}/Lock/BAM/${sample}_mapped.bam \
    -O ${batch_dir}/Lock/BAM/${sample}_dup_sort.bam \
    -M ${batch_dir}/Lock/BAM/${sample}_marked_dup_metrics.txt \
    --conf 'spark.executor.cores=${n_thread}' \
    --create-output-bam-index false \
    --create-output-bam-splitting-index false

  rm ${batch_dir}/Lock/BAM/${sample}_mapped.bam

  $picard ValidateSamFile I=${batch_dir}/Lock/BAM/${sample}_dup_sort.bam MODE=SUMMARY

  #########################################
  # Base quality score recalibration (BQSR):
  gatk BaseRecalibrator \
    -R $hg \
    --known-sites /home/projects/cu_10184/people/haikon/Reference/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /home/projects/cu_10184/people/haikon/Reference/GATK/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
    -I ${batch_dir}/Lock/BAM/${sample}_dup_sort.bam \
    -O ${batch_dir}/Lock/BAM/${sample}_recal_data.table

  gatk ApplyBQSR \
    -R $hg \
    -I ${batch_dir}/Lock/BAM/${sample}_dup_sort.bam \
    --bqsr-recal-file ${batch_dir}/Lock/BAM/${sample}_recal_data.table \
    -O ${batch_dir}/Lock/BAM/${sample}.bam

  rm ${batch_dir}/Lock/BAM/${sample}_dup_sort.bam
  mv ${batch_dir}/Lock/BAM/${sample}_marked_dup_metrics.txt ${batch_dir}/Lock/BAM/lock/
  mv ${batch_dir}/Lock/BAM/${sample}_recal_data.table ${batch_dir}/Lock/BAM/lock/
  #########################################
fi

####################################################################################################################
# Call SNV_InDel variants and annotate.
####################################################################################################################

# Annotation with Funcotator:

# Fill missing DP and AF in Funcotator maf output.


####################################################################################################################
####################################################################################################################
