####################################################################################################################
####################################################################################################################
# Filter SNV_InDel.
# Author: Haiying Kong
# Last Modified: 22 October 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:b:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    b) batch="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done

####################################################################################################################
####################################################################################################################
# Check if argument inputs from command are correct.
if [ -z "${dir_name}" ]
then
  echo "Error: Directory name is empty"
  exit 1
fi

if [ -z "$batch" ]
then
  echo "Error: Batch name is empty"
  exit 1
fi

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/Filter
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${batch_dir}/error/Filter
rm -rf ${error_dir}
mkdir -p ${error_dir}

# RedFlag directory:
rm -rf ${batch_dir}/RedFlag
mkdir -p ${batch_dir}/RedFlag

####################################################################################################################
####################################################################################################################
# Define directories.
####################################################################################################################
rm -rf ${batch_dir}/Result/SNV_InDel/Filtered
mkdir -p ${batch_dir}/Result/SNV_InDel/Filtered/Long
mkdir ${batch_dir}/Result/SNV_InDel/Filtered/Medium
mkdir ${batch_dir}/Result/SNV_InDel/Filtered/Short

####################################################################################################################
####################################################################################################################
# Get fastq file names.
cd ${batch_dir}/Result/SNV_InDel/AllVariants/Callers_Wide
maf_files=($(ls *.maf))

####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Get sample names.
samples=($(echo ${maf_files[@]%.maf} | tr ' ' '\n' | sort -u | tr '\n' ' '))

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
for sample in ${samples[@]}
do
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_Filter \
    -v dir_name=${dir_name},batch=${batch},sample=${sample},batch_dir=${batch_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/PDX/SNV_InDel/Filter_job.sh
done

####################################################################################################################
####################################################################################################################