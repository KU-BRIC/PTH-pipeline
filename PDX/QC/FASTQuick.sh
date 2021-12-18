####################################################################################################################
####################################################################################################################
# Quality control with QC.
# Author: Haiying Kong
# Last Modified: 22 October 2021
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
elif [ "$panel" = "panel1" ]
then 
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38
elif [ "$panel" = "panel2" ]
then 
  # Find target file for panel version 2.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2
elif [ "$panel" = "panel3" ]
then 
  # Find target file for panel version 3.
  target_name=Focused_myeloid_panel-All_target_segments_covered_by_probes-TE-93310852_hg38_v2_190722165759
else
  echo "Error: Please input panel version from command line as panel1, panel2 or panel3, or update BatchInfo.txt"
  exit 1
fi

if [ -z "${n_thread}" ]
then
  n_thread=8
fi

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
fq_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq/xengsort
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/FASTQuick
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${batch_dir}/error/FASTQuick
rm -rf ${error_dir}
mkdir -p ${error_dir}

# RedFlag directory:
rm -rf ${batch_dir}/RedFlag
mkdir -p ${batch_dir}/RedFlag

####################################################################################################################
####################################################################################################################
# QC:
####################################################################################################################
rm -rf ${batch_dir}/QC/FASTQuick
mkdir -p ${batch_dir}/QC/FASTQuick/Lock
mkdir -p ${batch_dir}/QC/FASTQuick/Result

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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_FASTQuick \
    -v n_thread=${n_thread},target_name=${target_name},dir_name=${dir_name},batch=${batch},sample=${sample},fq_dir=${fq_dir},batch_dir=${batch_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/PDX/QC/FASTQuick_job.sh
done

####################################################################################################################
####################################################################################################################
