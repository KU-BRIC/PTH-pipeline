####################################################################################################################
####################################################################################################################
# Perform sample quanlity control with FASTQuick.
# Author: Balthasar Schlotmann
# Last Modified: 10 June 2021
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
# Define directories.
qc_dir=/home/projects/cu_10184/projects/${dir_name}/QC

# Change to working directory.
temp_dir=/home/projects/cu_10184/projects/${dir_name}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

# log directory:
log_dir=${qc_dir}/log/FASTQuick
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${qc_dir}/error/FASTQuick
rm -rf ${error_dir}
mkdir -p ${error_dir}

# fastq directory:
fq_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq

# intermediate file directory:
lock_dir=${qc_dir}/Lock/FASTQuick
rm -rf ${lock_dir}
lock_dir=${lock_dir}/index/${batch}
mkdir -p ${lock_dir}

# result directory:
res_dir=${qc_dir}/Result/FASTQuick
rm -rf ${res_dir}
res_dir=${res_dir}/ByBatch/${batch}
mkdir -p ${res_dir}

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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_FASTQuick  \
    -v n_thread=${n_thread},target=${target},batch=${batch},sample=${sample},fq_dir=${fq_dir},lock_dir=${lock_dir},res_dir=${res_dir},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/QC/FASTQuick/FASTQuick_job.sh
done

####################################################################################################################
####################################################################################################################
