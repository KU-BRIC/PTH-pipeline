####################################################################################################################
####################################################################################################################
# Pre-process, call variants, annotate and filter variants.
# Author: Haiying Kong
# Last Modified: 18 July 2021
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
    echo "Error: sh /home/projects/cu_10184/projects/PTH/Code/Primary/Ensemble/Ensemble.sh -d PTH -b [batch_name] -p [panel_name] -t 8BatchInfo.txt does not have any information for this batch."
    exit 1
  fi
elif [ "$panel" = "panel1" ]
then
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38
elif [ "$panel" = "panel2" ]
then
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2
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
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}

# Create a new directory for the batch.
rm -rf ${batch_dir}
mkdir -p ${batch_dir}

# Change to working directory.
mkdir ${batch_dir}/temp
cd ${batch_dir}/temp

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
mkdir ${batch_dir}/log
log_dir=${batch_dir}/log/Ensemble
mkdir ${log_dir}

# error directory:
mkdir ${batch_dir}/error
error_dir=${batch_dir}/error/Ensemble
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
mkdir ${batch_dir}/Lock

####################################################################################################################
# Lock for BAM:
mkdir -p ${batch_dir}/Lock/BAM/lock

####################################################################################################################
# Lock for variant calling:
mkdir ${batch_dir}/Lock/SNV_InDel

####################################################################################################################
####################################################################################################################
# Define directories to save final results.
####################################################################################################################
mkdir ${batch_dir}/Result

####################################################################################################################
# Result for variant calling:
mkdir ${batch_dir}/Result/SNV_InDel

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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_Ensemble \
    -v n_thread=${n_thread},target_name=${target_name},dir_name=${dir_name},batch=${batch},sample=${sample},fq_dir=${fq_dir},batch_dir=${batch_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/Primary_1/Ensemble/Ensemble_job.sh
done

####################################################################################################################
####################################################################################################################
