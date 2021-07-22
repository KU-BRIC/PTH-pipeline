####################################################################################################################
####################################################################################################################
# Run MuTect2 on tumor samples, and annotate variants - job.
# Author: Haiying Kong
# Last Modified: 22 July 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=20:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
target="/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2.bed"
gnomad_DB="/home/projects/cu_10184/people/haikon/Reference/Funcotator/gnomAD_2.1_VCF_INFO_AF_Only/hg38/gnomad.genomes.r2.1.sites.liftoverToHg38.INFO_ANNOTATIONS_FIXED.vcf.gz"
pon=/home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon_1.vcf.gz
Funcotator_DB="/home/projects/cu_10184/people/haikon/Reference/Funcotator/funcotator_dataSources.v1.7.20200521s"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${batch_dir}/temp

####################################################################################################################
# Call SNV_InDels with MuTect2.
gatk Mutect2  \
  -R $hg -L ${target}  \
  --germline-resource ${gnomad_DB}  \
  --panel-of-normals ${pon}  \
  --base-quality-score-threshold 20 --min-base-quality-score 20  \
  -I ${batch_dir}/Lock/BAM/${sam}.bam  \
  -O ${batch_dir}/Lock/SNV_InDel/MuTect2_1/vcf/${sam}.vcf

####################################################################################################################
# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${batch_dir}/Lock/SNV_InDel/MuTect2_1/vcf/${sam}.vcf -O ${batch_dir}/Lock/SNV_InDel/MuTect2_1/maf/${sam}.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/MuTect2_1/Somatic/Fill_DP_AF_FuncotatorMAF.R ${batch_dir}/Lock/SNV_InDel/MuTect2_1 ${sam}

####################################################################################################################
####################################################################################################################
