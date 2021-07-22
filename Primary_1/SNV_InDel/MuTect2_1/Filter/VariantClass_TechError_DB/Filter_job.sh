####################################################################################################################
####################################################################################################################
# Filter SNV-InDels with filtering schemes - job.
# Author: Haiying Kong
# Last Modified: 21 July 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=10GB
#PBS -l walltime=1:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Filter.
####################################################################################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/Filter/VariantClass_TechError_DB/Filter_Variants.R ${batch_dir} ${sam}

####################################################################################################################
####################################################################################################################
