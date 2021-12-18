####################################################################################################################
####################################################################################################################
# Fingerprinting with Somalier to relate.
# Author: Haiying Kong
# Last Modified: 29 October 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

dir_name=/home/projects/cu_10184/projects/PTH/PAT_PDX
pedigree=${dir_name}/Meta/Pedigree.txt

####################################################################################################################
####################################################################################################################
# Create or clean the directories.
####################################################################################################################
# Lock directory:
lock_dir=${dir_name}/Lock/Somalier/Extracted

# Result directory:
res_dir=${dir_name}/Result/Somalier/
rm -rf ${res_dir}
mkdir -p ${res_dir}

# Working directory:
temp_dir=${dir_name}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

# log directory:
log_dir=${dir_name}/log/Somalier/Relate
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${dir_name}/error/Somalier/Relate
rm -rf ${error_dir}
mkdir -p ${error_dir}

# RedFlag directory:
rm -rf ${dir_name}/RedFlag
mkdir -p ${dir_name}/RedFlag

####################################################################################################################
####################################################################################################################
# Run Somalier relate.
qsub -o ${log_dir}/relate.log -e ${error_dir}/relate.error -N Somalier_Relate \
  -v lock_dir=${lock_dir},res_dir=${res_dir},temp_dir=${temp_dir} \
  /home/projects/cu_10184/projects/PTH/Code/PAT_PDX/Somalier/Somalier_Relate_job.sh

####################################################################################################################
####################################################################################################################
