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
    b) batch="$OPTARG";;
    p) panel="$OPTARG";;
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

if [ -z "$batch" ]
then
  echo "Error: Batch name is empty"
  exit 1
fi

if [ -z "$panel" ]
then
  # find target file for this batch from BatchInfo.txt with batch name.
  target_name=$(more /home/projects/cu_10184/projects/${dir_name}/Meta/BatchInfo.txt | awk -F '\t' -v batch="$batch" '( $1==batch ) {print $3}')
  if [ "${target_name}" = "" ]
  then
    echo "Error: BatchInfo.txt does not have any information for this batch."
    exit 1
  fi
  target=/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37/${target_name}.bed
elif [ "$panel" = "panel1" ]
then
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38
  target=/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37/${target_name}.bed
elif [ "$panel" = "panel2" ]
then
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2
  target=/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37/${target_name}.bed
else
  echo "Error: Please input panel version from command line as panel1 or panel2, or update BatchInfo.txt"
  exit 1
fi

if [ -z "${n_thread}" ]
then
  echo "By default, each job will use 8 cores."
  n_thread=8
fi 

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
fq_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq
batch_dir=/home/projects/cu_10184/projects/${dir_name}/QC/${batch}

# Create a new directory for the batch.
rm -rf ${batch_dir}
mkdir -p ${batch_dir}

# Change to working directory.
temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
mkdir ${batch_dir}/log
log_dir=${batch_dir}/log/FASTQuick_all
mkdir ${log_dir}

# error directory:
mkdir ${batch_dir}/error
error_dir=${batch_dir}/error/FASTQuick_all
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save QC results.
####################################################################################################################
QC_dir=${batch_dir}/FASTQuick_all
mkdir ${QC_dir}

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
# create file containing all sample fastq files
rm ${QC_dir}/fastqfile.tsv

for sample in ${samples[@]}
do
echo -e "${fq_dir}/${sample}_R1.fq.gz\t${fq_dir}/${sample}_R2.fq.gz" >> ${QC_dir}/fastqfile.tsv
done

####################################################################################################################
# start job

qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_FASTQuick  \
    -v n_thread=${n_thread},target=${target},batch=${batch},fq_dir=${fq_dir},QC_dir=${QC_dir},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/QC/FASTQuick/FASTQuick_AllSamples_job.sh

####################################################################################################################
####################################################################################################################
