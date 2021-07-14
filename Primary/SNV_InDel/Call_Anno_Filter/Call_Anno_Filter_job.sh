####################################################################################################################
####################################################################################################################
# Call variants, annotate and filter variants - job.
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
target_nochr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
target_chr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed

####################################################################################################################
# Software tools:
export PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin/:$PATH

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Call SNV_InDel variants and annotate.
####################################################################################################################
# VarDict:
################################################
# Set parameters.
allel_freq=0.005
map_qual=20
phred=20
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
vardict -G $hg -N ${sample} -b ${batch_dir}/Lock/BAM/${sample}.bam ${target_nochr} -c 1 -S 2 -E 3 -t  \
  -f ${allel_freq} -O ${map_qual} -q ${phred} -P ${var_pos} -o ${Qratio} -I ${indel_size}  \
  -M ${min_match} -L ${min_SV} -w ${ins_size} -W ${ins_SD} --nmfreq ${nmfreq} --mfreq ${mfreq}  \
  | teststrandbias.R \
  | var2vcf_valid.pl -N ${sample} -E -f ${allel_freq} > ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sample}.vcf

# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sample}.vcf -O ${batch_dir}/Lock/SNV_InDel/VarDict/maf/${sample}.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Fill_DP_AF_FuncotatorMAF_nonSNVer_job.R ${batch_dir}/Lock/SNV_InDel/VarDict ${sample} 'VarDict'

####################################################################################################################
# SNVer:
################################################
# Set parameters.
num_hap=2
het=0.001
map_qual=20
base_qual=20
str_bias=0.0001
fish_thresh=0.0001
p_thresh=0.05
min_read_strand=1
min_r_alt_ref=0.25

# Run SNVer.
conda activate
snver -r $hg -i ${batch_dir}/Lock/BAM/${sample}.bam -o ${batch_dir}/Lock/SNV_InDel/SNVer/vcf/${sample}  \
  -l ${target_chr} -n ${num_hap} -het ${het} -mq ${map_qual} -bq ${base_qual} -s ${str_bias}  \
  -f ${fish_thresh} -p ${p_thresh} -a ${min_read_strand} -b ${min_r_alt_ref}
conda deactivate

# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${batch_dir}/Lock/SNV_InDel/SNVer/vcf/${sample}.filter.vcf -O ${batch_dir}/Lock/SNV_InDel/SNVer/maf/${sample}.snv.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation
gatk Funcotator \
  -R $hg -V ${batch_dir}/Lock/SNV_InDel/SNVer/vcf/${sample}.indel.filter.vcf -O ${batch_dir}/Lock/SNV_InDel/SNVer/maf/${sample}.indel.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Fill_DP_AF_FuncotatorMAF_SNVer_job.R ${batch_dir}/Lock/SNV_InDel/SNVer ${sample} 'SNVer'

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
# Combine SNV_InDel variants from 3 callers and filter.
####################################################################################################################
# Combine variants from 3 callers.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Combine_Variants_AllCallers.R ${batch_dir}/Lock/SNV_InDel ${batch_dir}/Result/SNV_InDel/AllVariants ${batch} ${sample}

####################################################################################################################
# Filter.
####################################################################################################################
# Long:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "Long"

# Medium:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "Medium"

# Short:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "Short"

####################################################################################################################
####################################################################################################################
