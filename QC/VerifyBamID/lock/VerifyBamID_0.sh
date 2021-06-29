####################################################################################################################
####################################################################################################################
# Perform sequence and sample quanlity control.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 6 May 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:b:p:t:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    t) n_thread="$OPTARG";;
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

if [ -z "${n_thread}" ]
then
  echo "By default, each job will use 8 cores."
  n_thread=8
fi 

####################################################################################################################
####################################################################################################################

batch_dir=/home/projects/cu_10184/projects/PTH/QC/Lock/FingerPrinting/VerifyBamID/Result

# Create a new directory for the batch.
rm -rf ${batch_dir}
mkdir -p ${batch_dir}


while read f; do

batch=$(echo $f | cut -f2 -d' ') 
sample=$(echo $f | cut -f3 -d' ')
batch2=$(echo $f | cut -f4 -d' ')
sample2=$(echo $f | cut -f5 -d' ')

# Define directory for the batch.
bam_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}/Lock/BAM/

# Change to working directory.
temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
mkdir ${batch_dir}/log
log_dir=${batch_dir}/log/VerifyBamID
mkdir ${log_dir}

# error directory:
mkdir ${batch_dir}/error
error_dir=${batch_dir}/error/VerifyBamID
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save QC results.
####################################################################################################################
QC_dir=${batch_dir}/VerfyBamID
mkdir ${QC_dir}

####################################################################################################################
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################

qsub -o ${log_dir}/${sample}_vs_${sample2}.log -e ${error_dir}/${sample}_vs_${sample2}.error -N ${batch}_${sample}_vs_${sample2}_VerifyBamID  \
  -v n_thread=${n_thread},batch=${batch},sample=${sample},bam_dir=${bam_dir},QC_dir=${QC_dir},temp_dir=${temp_dir},sample2=${sample2},batch2=${batch2}  \
  /home/projects/cu_10184/projects/PTH/Code/QC/VerfyBamID/VerifyBamID_job.sh

done < /home/projects/cu_10184/projects/PTH/QC/Lock/FingerPrinting/VerifyBamID/PatientList_MultiSamples.txt
####################################################################################################################
####################################################################################################################
