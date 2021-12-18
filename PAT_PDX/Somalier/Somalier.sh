####################################################################################################################
####################################################################################################################
# Fingerprinting with Somalier.
# Author: Haiying Kong
# Last Modified: 24 October 2021
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
lock_dir=${dir_name}/Lock/Somalier
rm -rf ${lock_dir}
mkdir -p ${lock_dir}/Extracted

# Result directory:
res_dir=${dir_name}/Result/Somalier
rm -rf ${res_dir}
mkdir -p ${res_dir}

# Working directory:
temp_dir=${dir_name}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

# log directory:
log_dir=${dir_name}/log/Somalier
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${dir_name}/error/Somalier
rm -rf ${error_dir}
mkdir -p ${error_dir}

# RedFlag directory:
rm -rf ${dir_name}/RedFlag
mkdir -p ${dir_name}/RedFlag

####################################################################################################################
####################################################################################################################
# Read pedigree by line and submit job for fingerprinting.
sed 1d $pedigree | while IFS=$'\t' read -r -a one
do
  pat_vcf=${dir_name}/Lock/SNV_InDel/vcf/PAT/${one[2]}.vcf
  pdx_vcf=/home/projects/cu_10184/projects/PTH/BatchWork/${one[3]}/Lock/SNV_InDel/VarDict/vcf/${one[5]}.vcf
  res_name=${res_dir}/${one[5]}

  qsub -o ${log_dir}/${one[2]}_${one[5]}.log -e ${error_dir}/${one[2]}_${one[5]}.error -N ${one[2]}_${one[5]}_Somalier \
    -v pat_vcf=${pat_vcf},pdx_vcf=${pdx_vcf},dir_name=${dir_name},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/PAT_PDX/Somalier_job.sh
done

####################################################################################################################
####################################################################################################################
