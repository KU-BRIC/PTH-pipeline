####################################################################################################################
####################################################################################################################
# Fingerprinting with Somalier.
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
lock_dir=${dir_name}/Lock/Somalier
rm -rf ${lock_dir}
mkdir -p ${lock_dir}/Extracted

# Working directory:
temp_dir=${dir_name}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

# log directory:
log_dir=${dir_name}/log/Somalier/Extract
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${dir_name}/error/Somalier/Extract
rm -rf ${error_dir}
mkdir -p ${error_dir}

# RedFlag directory:
rm -rf ${dir_name}/RedFlag
mkdir -p ${dir_name}/RedFlag

####################################################################################################################
####################################################################################################################
# Read pedigree by line to create array with sample list.
i=0
bams=()
while IFS=$'\t' read -r -a one
do
  if [ $i -ne 0 ]
  then
    bams[${#bams[@]}]="/home/projects/cu_10184/projects/PTH/BatchWork/${one[0]}/Lock/BAM/${one[2]}.bam"
    bams[${#bams[@]}]="/home/projects/cu_10184/projects/PTH/BatchWork/${one[3]}/Lock/BAM/${one[5]}.bam"
    i=$((i+1))
  else
    i=$((i+1))
  fi
done < $pedigree

bam_files=$(echo "${bams[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')

for bam in ${bam_files}
do
  sam="$(basename -- $bam)"
  sam=${sam%.bam}
  qsub -o ${log_dir}/${sam}.log -e ${error_dir}/${sam}.error -N ${sam}_Somalier_Extract  \
    -v bam=${bam},sam=${sam},dir_name=${dir_name},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/PAT_PDX/Somalier/Somalier_Extract_job.sh
done

####################################################################################################################
####################################################################################################################
