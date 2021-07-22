####################################################################################################################
####################################################################################################################
# Filter SNV-InDels with filtering schemes.
# Author: Haiying Kong
# Last Modified: 21 July 2021
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
  echo "Error: Directory name is empty."
  exit 1
fi

if [ -z "$batch" ]
then
  echo "Error: Batch name is empty."
  exit 1
fi

scheme_name=VariantClass_TechError_DB

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/SNV_InDel/MuTect2/Filter/${scheme_name}
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${batch_dir}/error/SNV_InDel/MuTect2/Filter/${scheme_name}
rm -rf ${error_dir}
mkdir -p ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save results.
####################################################################################################################
rm -rf ${batch_dir}/Result/SNV_InDel/MuTect2
mkdir -p ${batch_dir}/Result/SNV_InDel/MuTect2
mkdir ${batch_dir}/Result/SNV_InDel/MuTect2/VariantClass
mkdir ${batch_dir}/Result/SNV_InDel/MuTect2/TechError
mkdir ${batch_dir}/Result/SNV_InDel/MuTect2/DB

####################################################################################################################
####################################################################################################################
# Get maf file names.
cd ${batch_dir}/Lock/SNV_InDel/MuTect2/maf
count=`ls -1 *.maf 2>/dev/null | wc -l`

if [[ $count = 0 ]]
then 
  exit
fi

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
for sam in ${samples[@]}
do
  qsub -o ${log_dir}/${sam}.log -e ${error_dir}/${sam}.error -N ${batch}_${sam}_SNV_InDel_Filter_${scheme_name}  \
    -v batch_dir=${batch_dir},sam=${sam},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2/Filter/${scheme_name}/Filter_job.sh
done

####################################################################################################################
####################################################################################################################
