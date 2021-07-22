####################################################################################################################
####################################################################################################################
# Run MuTect2 on tumor samples, and annotate variants - job.
# Author: Haiying Kong
# Last Modified: 21 July 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
gnomad_DB="/home/projects/cu_10184/people/haikon/Reference/Funcotator/gnomAD_2.1_VCF_INFO_AF_Only/hg38/gnomad.genomes.r2.1.sites.liftoverToHg38.INFO_ANNOTATIONS_FIXED.vcf.gz"
Funcotator_DB="/home/projects/cu_10184/people/haikon/Reference/Funcotator/funcotator_dataSources.v1.7.20200521s"
pon=/home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon_1.vcf.gz


batch_dir=/home/projects/cu_10184/projects/PTH/BatchWork_1/Primary_008
sam=Horizon-CMP016

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Call SNV_InDels with MuTect2.
gatk Mutect2  \
  -R $hg -I ${batch_dir}/Lock/BAM/${sam}.bam -O ${sam}.vcf  \
  -L /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2.bed  \
  --germline-resource ${gnomad_DB} --base-quality-score-threshold 20 --min-base-quality-score 20 --panel-of-normals /home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon_1.vcf.gz

####################################################################################################################
####################################################################################################################
