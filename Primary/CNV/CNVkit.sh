####################################################################################################################
####################################################################################################################
# Run CNVkit for all samples in the batch.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 22 June 2021
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
log_dir=${batch_dir}/log/CNVkit
rm -rf ${log_dir}
mkdir ${log_dir}

# error directory:
error_dir=${batch_dir}/error/CNVkit
rm -rf ${error_dir}
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
Lock_dir=${batch_dir}/Lock

####################################################################################################################
# Lock for BAM:
BAM_dir=${Lock_dir}/BAM

####################################################################################################################
# Lock for CNV:
Lock_CNVkit_dir=${Lock_dir}/CNV/CNVkit
rm -rf ${Lock_CNVkit_dir}
mkdir -p ${Lock_CNVkit_dir}

####################################################################################################################
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
for sample in ${BAM_dir}/*.bam
do
  sample=$(basename $sample)
  sample=${sample/.bam/}
  echo $sample
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_CNVkit \
    -v batch=${batch},sample=${sample},BAM_dir=${BAM_dir},Lock_CNVkit_dir=${Lock_CNVkit_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Primary/CNV/CNVkit_job.sh
done

####################################################################################################################
####################################################################################################################
