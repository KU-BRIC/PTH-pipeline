####################################################################################################################
####################################################################################################################
# Filter SNV_InDel - job.
# Author: Haiying Kong
# Last Modified: 22 October 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=10GB
#PBS -l walltime=200:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Filter.
####################################################################################################################
# Long:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "PTH" "Long"

if [ $? -ne 0 ];  then echo "Failed at filtering with Long Scheme." >> ${batch_dir}/RedFlag/${sample}.txt;  fi

# Medium:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "PTH" "Medium"

if [ $? -ne 0 ];  then echo "Failed at filtering with Medium Scheme." >> ${batch_dir}/RedFlag/${sample}.txt;  fi

# Short:
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/SNV_InDel/Filter/Filter_Variants.R ${batch_dir}/Result/SNV_InDel ${sample} "PTH" "Short"

if [ $? -ne 0 ];  then echo "Failed at filtering with Short Scheme." >> ${batch_dir}/RedFlag/${sample}.txt;  fi

####################################################################################################################
####################################################################################################################
