####################################################################################################################
####################################################################################################################
# Identify ITD with getITD - job.
# Author: Haiying Kong
# Last Modified: 23 April 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=50:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
reference="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon.txt"
anno="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv"

# Software tools:
getITD="python3 /home/projects/cu_10184/people/haikon/Software/getITD/getitd.py"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run getITD:
cd ${Lock_getITD_dir}
conda activate getITD
$getITD -reference $reference -anno $anno -infer_sense_from_alignment True -plot_coverage True -require_indel_free_primers False -min_read_length 70 -min_bqs 20 -min_read_copies 1 -filter_ins_vaf 0.001 ${sample} ${fq_dir}/${sample}_R1.fq.gz ${fq_dir}/${sample}_R2.fq.gz
conda deactivate

####################################################################################################################
# Clean the output folder.
mv ${sample}_getitd ${sample}

####################################################################################################################
####################################################################################################################
