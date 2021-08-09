####################################################################################################################
####################################################################################################################
# Pre-process, call variants, annotate and filter variants - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 8 August 2021
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

# CNVkit PoN:
cnvkit_pon="/home/projects/cu_10184/projects/PTH/Reference/CNVkit/OutputReference.cnn"

# ITD:
target_flt3="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/FLT3.bed"
target_itd=/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/${target_name}.bed
getITD_reference="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon.txt"
getITD_anno="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv"

# FASTQuick:
fastquick_reference_dir="/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37"
fastquick_resource_dir="/home/projects/cu_10184/people/haikon/Software/FASTQuick/resource"
target_fastquick=/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37/${target_name}.bed

####################################################################################################################
# Software tools:
picard="java -jar /home/projects/cu_10184/people/haikon/Software/picard.jar"
ScanITD="/home/projects/cu_10184/people/haikon/Software/ScanITD/ScanITD.py"
getITD="python3 /home/projects/cu_10184/people/haikon/Software/getITD/getitd.py"
FASTQuick="/home/projects/cu_10184/people/haikon/Software/FASTQuick/bin/FASTQuick.sh"

export PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin/:$PATH

####################################################################################################################
# Set parameters.
itd_scheme="Scheme_2"
itd_thresh_n_alt=1

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
  | var2vcf_valid.pl -N ${sample} -E -f ${allel_freq} > ${batch_dir}/Lock/SNV_InDel/VarDict/vcf_0/${sample}.vcf

# Filter for PASS.
awk -F '\t' '{if($1~/^#/ || $7=="PASS") print $0}' ${batch_dir}/Lock/SNV_InDel/VarDict/vcf_0/${sample}.vcf > ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sample}.vcf

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
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "PTH" "Long"

# Medium:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "PTH" "Medium"

# Short:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "PTH" "Short"

####################################################################################################################
####################################################################################################################
# Call CNVs.
####################################################################################################################
# CNACS:
conda activate CNACS

export JAVAPATH=/home/projects/cu_10184/people/haikon/Software/anaconda3/bin
export PICARD_PATH=/home/projects/cu_10184/people/haikon/Software/picard.jar
export SAMTOOLS_PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin
export PERL_PATH=/usr/bin/perl
export BEDTOOLS_PATH=/home/projects/cu_10184/people/haikon/Software/anaconda3/bin
export R_PATH=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/R
export R_LIBS_PATH=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/library
export R_LIBS=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/library

rm -rf ${batch_dir}/Lock/CNV/CNACS/${sample}
mkdir -p ${batch_dir}/Lock/CNV/CNACS/${sample}

toil_cnacs run \
    ${batch_dir}/Lock/CNV/CNACS/${sample}/jobstore/ \
    --stats \
    --db_dir /home/projects/cu_10184/projects/PTH/Reference/CNACS/db \
    --outdir ${batch_dir}/Lock/CNV/CNACS \
    --pool_dir /home/projects/cu_10184/projects/PTH/Reference/CNACS/PoN \
    --fasta $hg \
    --samp ${batch_dir}/Lock/BAM/${sample}.bam

rm -r ${batch_dir}/Lock/CNV/CNACS/${sample}/tmp
rm -r ${batch_dir}/Lock/CNV/CNACS/${sample}/jobstore

conda deactivate

####################################################################################################################
# CNVkit:
conda activate CNVkit

cnvkit.py batch ${batch_dir}/Lock/BAM/${sample}.bam \
    --reference ${cnvkit_pon} \
    --output-dir ${batch_dir}/Lock/CNV/CNVkit

echo "cnvkit.py scatter ${batch_dir}/Lock/CNV/CNVkit/${sample}.cnr -s ${batch_dir}/Lock/CNV/CNVkit/${sample}.cns -o ${batch_dir}/Lock/CNV/CNVkit/${sample}.scatter.pdf"

cnvkit.py scatter ${batch_dir}/Lock/CNV/CNVkit/${sample}.cnr -s ${batch_dir}/Lock/CNV/CNVkit/${sample}.cns -o ${batch_dir}/Lock/CNV/CNVkit/${sample}.scatter.pdf
cnvkit.py diagram ${batch_dir}/Lock/CNV/CNVkit/${sample}.cnr -s ${batch_dir}/Lock/CNV/CNVkit/${sample}.cns -o ${batch_dir}/Lock/CNV/CNVkit/${sample}.diagram.pdf

conda deactivate

####################################################################################################################
####################################################################################################################
# Get depth of coverage profile.
####################################################################################################################
mkdir -p ${batch_dir}/Lock/DepthOfCoverage/FreqTable
mkdir -p ${batch_dir}/Lock/DepthOfCoverage/DensityPlot

Rscript /home/projects/cu_10184/projects/PTH/Code/Source/DOC/HelloRanges_one_sample_CoverageFreqPlot.R ${batch_dir}/Lock/BAM ${sample}.bam ${target_nopad} ${batch_dir}/Lock/DepthOfCoverage

####################################################################################################################
####################################################################################################################
# ITD identification on FLT3.
####################################################################################################################
# VarDict:
################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${itd_scheme}/VarDict_FilterClean.R ${batch} ${sample} ${batch_dir}/Lock/SNV_InDel/VarDict ${batch_dir}/Lock/ITD/VarDict

####################################################################################################################
# Pindel:
################################################
# Run Pindel.
conda activate
echo -e ${batch_dir}/Lock/BAM/${sample}.bam"\t"250"\t"${sample} > ${batch_dir}/Lock/ITD/Pindel/${sample}_config.txt
pindel -f $hg -t ${target_itd} -c chr13:28033736-28034557  \
       -i ${batch_dir}/Lock/ITD/Pindel/${sample}_config.txt -o ${batch_dir}/Lock/ITD/Pindel/${sample}  \
       -T $n_thread -x 3 -u 0.04
rm ${batch_dir}/Lock/ITD/Pindel/${sample}_config.txt
pindel2vcf -p ${batch_dir}/Lock/ITD/Pindel/${sample}_TD -r $hg -R HG38 -d 20201224 -v ${batch_dir}/Lock/ITD/Pindel/${sample}.vcf
conda deactivate

# Clean the output vcf.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${itd_scheme}/Pindel_Clean.R ${batch} ${sample} ${batch_dir}/Lock/ITD/Pindel

####################################################################################################################
# ScanITD:
################################################
# Run ScanITD.
conda activate ScanITD
$ScanITD -r $hg -t ${target_itd} -i ${batch_dir}/Lock/BAM/${sample}.bam -o ${batch_dir}/Lock/ITD/ScanITD/${sample} -m 20 -c 1 -d 100 -f 0.0001 -l 5 -n 3
conda deactivate

#  -m MAPQ, --mapq MAPQ  minimal MAPQ in BAM for calling ITD (default: 15)
#  -c AO, --ao AO        minimal observation count for ITD (default: 4)
#  -d DP, --depth DP     minimal depth to call ITD (default: 10)
#  -f VAF, --vaf VAF     minimal variant allele frequency (default: 0.1)
#  -l ITD_LEN, --len ITD_LEN
#                        minimal ITD length to report (default: 10)
#  -n MISMATCH           maximum allowed mismatch bases of pairwise local alignment (default: 3)

# Clean the output tsv.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${itd_scheme}/ScanITD_Clean.R ${batch} ${sample} ${batch_dir}/Lock/ITD/ScanITD

####################################################################################################################
# getITD:
################################################
# Filter bam file with target bed file and conver it to fastq.
mkdir ${batch_dir}/Lock/ITD/getITD/${sample}
cd ${batch_dir}/Lock/ITD/getITD/${sample}
bedtools intersect -abam ${batch_dir}/Lock/BAM/${sample}.bam -b ${target_flt3} -u > FLT3.bam
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
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${itd_scheme}/getITD_Clean.R ${batch} ${sample} ${batch_dir}/Lock/ITD/getITD

####################################################################################################################
# Plot zoom out IGV for whole FLT3 region.
################################################
rm -rf ${batch_dir}/Lock/BAM/${sample}.bam.bai
ln -s ${batch_dir}/Lock/BAM/${sample}.bai ${batch_dir}/Lock/BAM/${sample}.bam.bai

cat > ${batch_dir}/Lock/ITD/IGV/${sample}_batch.bat << EOF
new 
genome hg38 
snapshotDirectory ${batch_dir}/Lock/ITD/IGV 
load ${batch_dir}/Lock/BAM/${sample}.bam
maxPanelHeight 20000
preference SAM.SHOW_SOFT_CLIPPED TRUE
goto chr13:28033736-28034557
snapshot ${sample}.png 
exit
EOF

/home/projects/cu_10184/people/haikon/Software/IGV_Linux_2.9.0/igv_auto.sh -b ${batch_dir}/Lock/ITD/IGV/${sample}_batch.bat
rm ${batch_dir}/Lock/ITD/IGV/${sample}_batch.bat

####################################################################################################################
# Combine the outcome from different tools to come up with final ITD list.
################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${itd_scheme}/Combine_Filter.R ${batch} ${sample} ${batch_dir}/Lock/ITD ${batch_dir}/Result/ITD ${itd_thresh_n_alt}

####################################################################################################################
# Plot zoom in IGV if any ITD is identified for this sample.
################################################
bat_files=${batch_dir}/Result/ITD/IGV/${sample}*.bat

shopt -s nullglob
for batfile in ${batch_dir}/Result/ITD/IGV/${sample}*.bat
do
  /home/projects/cu_10184/people/haikon/Software/IGV_Linux_2.9.0/igv_auto.sh -b ${batfile}
  cat $batfile
  rm ${batfile}
done 

# Remove .bam.bai and .bat files.
rm ${batch_dir}/Lock/BAM/${sample}.bam.bai

####################################################################################################################
####################################################################################################################
# Quality Control.
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
