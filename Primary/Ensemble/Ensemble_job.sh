####################################################################################################################
####################################################################################################################
# Pre-process, call variants, annotate and filter variants - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 19 May 2021
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
target_flt3="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/FLT3.bed"
getITD_reference="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon.txt"
getITD_anno="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv"

# Software tools:
picard="java -jar /home/projects/cu_10184/people/haikon/Software/picard.jar"
getITD="python3 /home/projects/cu_10184/people/haikon/Software/getITD/getitd.py"
       
# Set parameters.
thresh_n_alt=1

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Preprocessing with BWA, Picard and GATK.
####################################################################################################################
# Map to hg38:
bwa mem \
  $hg \
  ${fq_dir}/${sample}_R1.fq.gz ${fq_dir}/${sample}_R2.fq.gz \
  -t ${n_thread} \
  > ${BAM_dir}/${sample}_bwa.sam

$picard FastqToSam \
  F1=${fq_dir}/${sample}_R1.fq.gz \
  F2=${fq_dir}/${sample}_R2.fq.gz \
  O=${BAM_dir}/${sample}_unmapped.sam \
  PL=ILLUMINA \
  SM=${sample}

$picard MergeBamAlignment \
  R=$hg \
  ALIGNED=${BAM_dir}/${sample}_bwa.sam \
  UNMAPPED=${BAM_dir}/${sample}_unmapped.sam \
  O=${BAM_dir}/${sample}_mapped.bam

rm ${BAM_dir}/${sample}_bwa.sam
rm ${BAM_dir}/${sample}_unmapped.sam

#########################################
# Mark duplicates and sort:
gatk MarkDuplicatesSpark \
  -I ${BAM_dir}/${sample}_mapped.bam \
  -O ${BAM_dir}/${sample}_dup_sort.bam \
  -M ${BAM_dir}/${sample}_marked_dup_metrics.txt \
  --conf 'spark.executor.cores=${n_thread}' \
  --create-output-bam-index false \
  --create-output-bam-splitting-index false

rm ${BAM_dir}/${sample}_mapped.bam

$picard ValidateSamFile I=${BAM_dir}/${sample}_dup_sort.bam MODE=SUMMARY

#########################################
# Base quality score recalibration (BQSR):
gatk BaseRecalibrator \
  -R $hg \
  --known-sites /home/projects/cu_10184/people/haikon/Reference/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites /home/projects/cu_10184/people/haikon/Reference/GATK/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  -I ${BAM_dir}/${sample}_dup_sort.bam \
  -O ${BAM_dir}/${sample}_recal_data.table

gatk ApplyBQSR \
  -R $hg \
  -I ${BAM_dir}/${sample}_dup_sort.bam \
  --bqsr-recal-file ${BAM_dir}/${sample}_recal_data.table \
  -O ${BAM_dir}/${sample}.bam
 
rm ${BAM_dir}/${sample}_dup_sort.bam
mv ${BAM_dir}/${sample}_marked_dup_metrics.txt ${BAM_lock_dir}/
mv ${BAM_dir}/${sample}_recal_data.table ${BAM_lock_dir}/

####################################################################################################################
# Call SNV_InDel variants and annotate.
####################################################################################################################
# VarDict:
################################################
# Set parameters.
allel_freq=0.01
map_qual=20
phred=22.5
var_pos=5
Qratio=1.5
indel_size=50
min_match=25
min_SV=500
ins_size=250
ins_SD=100
nmfreq=0.1
mfreq=0.25

# Run VarDict:
vardict -G $hg -N ${sample} -b ${BAM_dir}/${sample}.bam $target_nochr -c 1 -S 2 -E 3 -t  \
  -f ${allel_freq} -O ${map_qual} -q ${phred} -P ${var_pos} -o ${Qratio} -I ${indel_size}  \
  -M ${min_match} -L ${min_SV} -w ${ins_size} -W ${ins_SD} --nmfreq ${nmfreq} --mfreq ${mfreq}  \
  | teststrandbias.R \
  | var2vcf_valid.pl -N ${sample} -E -f ${allel_freq} > ${Lock_VarDict_dir}/vcf/${sample}.vcf

# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${Lock_VarDict_dir}/vcf/${sample}.vcf -O ${Lock_VarDict_dir}/maf/${sample}.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Fill_DP_AF_FuncotatorMAF_nonSNVer_job.R ${Lock_VarDict_dir} ${sample} 'VarDict'

####################################################################################################################
# SNVer:
################################################
# Set parameters.
num_hap=2
het=0.001
map_qual=20
base_qual=17
str_bias=0.0001
fish_thresh=0.0001
p_thresh=0.05
min_read_strand=1
min_r_alt_ref=0.25

# Run SNVer.
conda activate
snver -r $hg -i ${BAM_dir}/${sample}.bam -o ${Lock_SNVer_dir}/vcf/${sample}  \
  -l ${target_chr} -n ${num_hap} -het ${het} -mq ${map_qual} -bq ${base_qual} -s ${str_bias}  \
  -f ${fish_thresh} -p ${p_thresh} -a ${min_read_strand} -b ${min_r_alt_ref}
conda deactivate

# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${Lock_SNVer_dir}/vcf/${sample}.filter.vcf -O ${Lock_SNVer_dir}/maf/${sample}.snv.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation
gatk Funcotator \
  -R $hg -V ${Lock_SNVer_dir}/vcf/${sample}.indel.filter.vcf -O ${Lock_SNVer_dir}/maf/${sample}.indel.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Fill_DP_AF_FuncotatorMAF_SNVer_job.R ${Lock_SNVer_dir} ${sample} 'SNVer'

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
# Combine SNV_InDel variants from 3 callers and filter.
####################################################################################################################
# Combine variants from 3 callers.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Combine_Variants_AllCallers.R ${Lock_SNV_InDel_dir} ${Result_SNV_InDel_dir}/AllVariants ${batch} ${sample}

####################################################################################################################
# Filter.
####################################################################################################################
# Long:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants_Long.R ${Result_SNV_InDel_dir} ${sample}

# Medium:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants_Medium.R ${Result_SNV_InDel_dir} ${sample}

# Short:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants_Short.R ${Result_SNV_InDel_dir} ${sample}

####################################################################################################################
####################################################################################################################
# Call CNVs and annotate.
####################################################################################################################
conda activate CNACS

export JAVAPATH=/home/projects/cu_10184/people/haikon/Software/anaconda3/bin
export PICARD_PATH=/home/projects/cu_10184/people/haikon/Software/picard.jar
export SAMTOOLS_PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin
export PERL_PATH=/usr/bin/perl
export BEDTOOLS_PATH=/home/projects/cu_10184/people/haikon/Software/anaconda3/bin
export R_PATH=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/R
export R_LIBS_PATH=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/library
export R_LIBS=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/library

rm -rf ${Lock_CNACS_dir}/${sample}
mkdir -p ${Lock_CNACS_dir}/${sample}

toil_cnacs run \
    ${Lock_CNACS_dir}/${sample}/jobstore/ \
    --stats \
    --db_dir /home/projects/cu_10184/projects/PTH/Reference/CNACS/db \
    --outdir ${Lock_CNACS_dir} \
    --pool_dir /home/projects/cu_10184/projects/PTH/Reference/CNACS/PoN \
    --fasta $hg \
    --samp ${BAM_dir}/${sample}.bam

rm -r ${Lock_CNACS_dir}/${sample}/tmp
rm -r ${Lock_CNACS_dir}/${sample}/jobstore

conda deactivate

####################################################################################################################
####################################################################################################################
# Get depth of coverage profile.
####################################################################################################################
mkdir -p ${Lock_DOC_dir}/FreqTable
mkdir -p ${Lock_DOC_dir}/DensityPlot

Rscript /home/projects/cu_10184/projects/PTH/Code/Source/DOC/HelloRanges_one_sample_CoverageFreqPlot.R ${BAM_dir} ${sample}.bam ${target_nopad} ${Lock_DOC_dir}

####################################################################################################################
####################################################################################################################
# ITD identification on FLT3.
####################################################################################################################
####################################################################################################################
# VarDict:
####################################################################################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/VarDict_FilterClean.R ${batch} ${sample} ${Lock_VarDict_dir} ${Lock_ITD_dir}/VarDict

####################################################################################################################
# Pindel:
####################################################################################################################
# Run Pindel.
conda activate
echo -e ${BAM_dir}/${sample}.bam"\t"250"\t"${sample} > ${Lock_ITD_dir}/Pindel/${sample}_config.txt
pindel -f $hg -t ${target_itd} -c chr13:28033736-28034557  \
       -i ${Lock_ITD_dir}/Pindel/${sample}_config.txt -o ${Lock_ITD_dir}/Pindel/${sample}  \
       -T $n_thread -x 3 -u 0.04
rm ${Lock_ITD_dir}/Pindel/${sample}_config.txt
pindel2vcf -p ${Lock_ITD_dir}/Pindel/${sample}_TD -r $hg -R HG38 -d 20201224 -v ${Lock_ITD_dir}/Pindel/${sample}.vcf
conda deactivate

# Clean the output vcf.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/Pindel_Clean.R ${batch} ${sample} ${Lock_ITD_dir}/Pindel

####################################################################################################################
# getITD:
####################################################################################################################
# Filter bam file with target bed file and conver it to fastq.
mkdir ${Lock_ITD_dir}/getITD/${sample}
cd ${Lock_ITD_dir}/getITD/${sample}
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
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/getITD_Clean.R ${batch} ${sample} ${Lock_ITD_dir}/getITD

####################################################################################################################
# Combine and the results to come up with final ITD list.
####################################################################################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/CombineFilter_Pindel_getITD_VarDictAnno.R ${batch} ${sample} ${Lock_ITD_dir} ${Result_ITD_dir} ${thresh_n_alt}

####################################################################################################################
####################################################################################################################
