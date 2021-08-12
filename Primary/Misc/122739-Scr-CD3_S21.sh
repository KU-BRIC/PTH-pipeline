####################################################################################################################
####################################################################################################################
# Re-run from Funcotator on VarDict output, combines results from callers, filter, CNV, ITD.
# Author: Haiying Kong
# Last Modified: 9 August 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

dir_name=PTH_test
batch=Kirsten_groenbaek_small_panel_data
batch_dir=/home/projects/cu_10184/projects/PTH_test/BatchWork/Kirsten_groenbaek_small_panel_data
temp_dir=/home/projects/cu_10184/projects/PTH_test/BatchWork/Kirsten_groenbaek_small_panel_data/temp
sam=122739-Scr-CD3_S21

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
Funcotator_DB="/home/projects/cu_10184/people/haikon/Reference/Funcotator/funcotator_dataSources.v1.7.20200521s"

# Panel target:
target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38
target_nochr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
target_chr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed
target_nopad=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Original/${target_name}.bed

####################################################################################################################
# Software tools:
export PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin/:$PATH

####################################################################################################################
####################################################################################################################
# Delete rows of problem.
mv ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sam}.vcf ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sam}_0.vcf
rm -f ${batch_dir}/Lock/SNV_InDel/VarDict/maf/${sam}.maf
awk -F '\t' '{if (!($1=="chr17" && ($2==76736873 || $2==76736971))) print $0}' ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sam}_0.vcf > ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sam}.vcf

####################################################################################################################
####################################################################################################################
# Script for th job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Call SNV_InDel variants and annotate.
####################################################################################################################
# VarDict:
################################################
# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${batch_dir}/Lock/SNV_InDel/VarDict/vcf/${sam}.vcf -O ${batch_dir}/Lock/SNV_InDel/VarDict/maf/${sam}.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Fill_DP_AF_FuncotatorMAF_nonSNVer_job.R ${batch_dir}/Lock/SNV_InDel/VarDict ${sam} 'VarDict'

####################################################################################################################
# Combine SNV_InDel variants from 3 callers and filter.
####################################################################################################################
# Combine variants from 3 callers.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Combine_Variants_AllCallers.R ${batch_dir}/Lock/SNV_InDel ${batch_dir}/Result/SNV_InDel/AllVariants ${batch} ${sam}

####################################################################################################################
# Filter.
####################################################################################################################
# Long:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sam} "PTH" "Long"

# Medium:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sam} "PTH" "Medium"

# Short:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sam} "PTH" "Short"

####################################################################################################################
####################################################################################################################
# ITD identification on FLT3.
####################################################################################################################
# VarDict:
################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${itd_scheme}/VarDict_FilterClean.R ${batch} ${sam} ${batch_dir}/Lock/SNV_InDel/VarDict ${batch_dir}/Lock/ITD/VarDict

####################################################################################################################
# Combine the outcome from different tools to come up with final ITD list.
################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${itd_scheme}/Combine_Filter.R ${batch} ${sam} ${batch_dir}/Lock/ITD ${batch_dir}/Result/ITD ${itd_thresh_n_alt}

####################################################################################################################
# Plot zoom in IGV if any ITD is identified for this sample.
################################################
rm -rf ${batch_dir}/Lock/BAM/${sam}.bam.bai
ln -s ${batch_dir}/Lock/BAM/${sam}.bai ${batch_dir}/Lock/BAM/${sam}.bam.bai

bat_files=${batch_dir}/Result/ITD/IGV/${sam}*.bat

shopt -s nullglob
for batfile in ${batch_dir}/Result/ITD/IGV/${sam}*.bat
do
  /home/projects/cu_10184/people/haikon/Software/IGV_Linux_2.9.0/igv_auto.sh -b ${batfile}
  cat $batfile
  rm ${batfile}
done 

# Remove .bam.bai and .bat files.
rm ${batch_dir}/Lock/BAM/${sam}.bam.bai

####################################################################################################################
####################################################################################################################
