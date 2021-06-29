####################################################################################################################
####################################################################################################################
# Perform sample quanlity control with FASTQuick.
# Author: Balthasar Schlotmann
# Last Modified: 24 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:b:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done
n_thread=${n_thread#0}

####################################################################################################################
####################################################################################################################
# Check if argument inputs from command are correct.
if [ -z "${dir_name}" ]
then
  echo "Error: Directory name is empty"
  exit 1
fi

####################################################################################################################
####################################################################################################################
# Define directories.
qc_dir=/home/projects/cu_10184/projects/${dir_name}/QC

# Change to working directory.
temp_dir=/home/projects/cu_10184/projects/${dir_name}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

# log directory:
log_dir=${qc_dir}/log/Somalier
mkdir -p ${log_dir}

# error directory:
error_dir=${qc_dir}/error/Somalier
mkdir -p ${error_dir}

# intermediate file directory:
lock_dir=${qc_dir}/Lock/Somalier/
mkdir -p ${lock_dir}

# result directory:
res_dir=${qc_dir}/Result/Somalier
rm -rf ${res_dir}
mkdir -p ${res_dir}

####################################################################################################################
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

qsub -o ${log_dir}/Relate.log -e ${error_dir}/Relate.error -N Relate_Somalier  \
  -v lock_dir=${lock_dir},res_dir=${res_dir}  \
  /home/projects/cu_10184/projects/PTH/Code/QC/Somalier/Somalier_relate_job.sh

####################################################################################################################
####################################################################################################################
