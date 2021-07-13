####################################################################################################################
####################################################################################################################
# Filter SNV-InDels with filtering schemes.
# Author: Haiying Kong
# Last Modified: 13 July 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:b:r:f:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    b) batch="$OPTARG";;
    r) scheme_dir="$OPTARG";;
    f) scheme_name="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done

####################################################################################################################
####################################################################################################################
# Check if argument inputs from command are correct.
if [ -z "${dir_name}" ]
then
  dir_name=PTH
  echo "Error: Directory name is empty."
  exit 1
fi

if [ -z "$batch" ]
then
  echo "Error: Batch name is empty."
  exit 1
fi

if [ -z "${scheme_dir}" ]
then
  scheme_name=PTH
  echo "By default, the filtering scheme is from the project directory PTH."
fi

if [ -z "${scheme_name}" ]
then
  scheme_name=NewScheme
  echo "By default, filtering scheme name is NewScheme."
fi

if [ ! -d "/home/projects/cu_10184/projects/${scheme_dir}/Reference/Filtering/${scheme_name}" ]
then
  echo "The reference files for the filtering scheme do not exist."
  exit 1
fi

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/SNV_InDel/Filter/${scheme_name}
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${batch_dir}/error/SNV_InDel/Filter/${scheme_name}
rm -rf ${error_dir}
mkdir -p ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save final results.
####################################################################################################################
Result_dir=${batch_dir}/Result/SNV_InDel
rm -rf ${Result_dir}/Filtered/${scheme_name}
mkdir -p ${Result_dir}/Filtered/${scheme_name}

####################################################################################################################
####################################################################################################################
# Get maf file names.
cd ${Result_dir}/AllVariants/Callers_Wide
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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_SNV_InDel_Filter_${scheme_name} \
    -v Result_dir=${Result_dir},sample=${sample},scheme_dir=${scheme_dir},scheme_name=${scheme_name},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Filter_${scheme_name}/Filter_job.sh
done

####################################################################################################################
####################################################################################################################
