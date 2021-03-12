####################################################################################################################
####################################################################################################################
# Classify xenofraft sample reads to human and mouse.
# Author: Haiying Kong
# Last Modified: 12 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

source /home/people/haikon/.bashrc

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:b:t:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    b) batch="$OPTARG";;
    t) n_thread="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&3;;
  esac
done
n_thread=${n_thread#0}

####################################################################################################################
####################################################################################################################
# Check if arguments are input from command.
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

if [ -z "${n_thread}" ]
then
  n_thread=8
fi

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
in_dir=/home/projects/cu_10184/projects/${dir_name}/PDXseqData/${batch}/fastq
out_dir=/home/projects/cu_10184/projects/${dir_name}/PDXseqData/${batch}/xengsort
human_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq/xengsort
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

####################################################################################################################
# Create new directories to save classified reads.
rm -rf ${out_dir}
mkdir -p ${out_dir}
mkdir ${out_dir}/Human
mkdir ${out_dir}/Mouse
mkdir ${out_dir}/Both
mkdir ${out_dir}/Neither
mkdir ${out_dir}/Ambiguous
mkdir ${out_dir}/ReadCount

# Create new directory to soft link human reads.
rm -rf ${human_dir}
mkdir -p ${human_dir}

# Create new directories to save log, error, and possible temp files.
mkdir -p ${batch_dir}

# log directory:
mkdir -p ${batch_dir}/log
log_dir=${batch_dir}/log/xengsort
rm -rf ${log_dir}
mkdir ${log_dir}

# error directory:
mkdir -p ${batch_dir}/error
error_dir=${batch_dir}/error/xengsort
rm -rf ${error_dir}
mkdir ${error_dir}

temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}

# Change to work directory.
cd ${temp_dir}

####################################################################################################################
# Get fastq file names.
cd ${in_dir}
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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_xengsort \
    -v n_thread=${n_thread},sample=${sample},in_dir=${in_dir},out_dir=${out_dir},human_dir=${human_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/PDX/Classify_Reads/xengsort_job.sh
done

####################################################################################################################
####################################################################################################################
