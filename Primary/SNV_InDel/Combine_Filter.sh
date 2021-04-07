####################################################################################################################
####################################################################################################################
# Combine variants called by 3 callers and filter them.
# Author: Haiying Kong
# Last Modified: 7 April 2021
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
fq_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/SNV_InDel/Combine_Filter
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${batch_dir}/error/SNV_InDel/Combine_Filter
rm -rf ${error_dir}
mkdir -p ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
Lock_dir=${batch_dir}/Lock

####################################################################################################################
# Lock for variant calling:
Lock_SNV_InDel_dir=${Lock_dir}/SNV_InDel

####################################################################################################################
####################################################################################################################
# Define directories to save final results.
####################################################################################################################
Result_dir=${batch_dir}/Result

####################################################################################################################
# Result for variant calling:
Result_SNV_InDel_dir=${Result_dir}/SNV_InDel
rm -rf ${Result_SNV_InDel_dir}
mkdir ${Result_SNV_InDel_dir}
mkdir ${Result_SNV_InDel_dir}/AllVariants
mkdir ${Result_SNV_InDel_dir}/AllVariants/Callers_Long
mkdir ${Result_SNV_InDel_dir}/AllVariants/Callers_Wide
mkdir ${Result_SNV_InDel_dir}/Filtered
mkdir ${Result_SNV_InDel_dir}/Filtered/Long
mkdir ${Result_SNV_InDel_dir}/Filtered/Medium
mkdir ${Result_SNV_InDel_dir}/Filtered/Short

####################################################################################################################
####################################################################################################################
# Get fastq file names.
cd ${fq_dir}
fq_files=($(ls *.fq.gz))

####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Get sample names.
samples=($(echo ${fq_files[@]%_R*.fq.gz} | tr ' ' '\n' | sort -u | tr '\n' ' '))

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
for sample in ${samples[@]}
do
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_SNV_InDel_Combine_Filter \
    -v batch=${batch},sample=${sample},Lock_SNV_InDel_dir=${Lock_SNV_InDel_dir},Result_SNV_InDel_dir=${Result_SNV_InDel_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Combine_Filter_job.sh
done

####################################################################################################################
####################################################################################################################
