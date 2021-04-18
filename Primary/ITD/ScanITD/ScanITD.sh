####################################################################################################################
####################################################################################################################
# Identify ITD with ScanITD.
# Author: Haiying Kong
# Last Modified: 18 April 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:b:p:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    b) batch="$OPTARG";;
    p) panel="$OPTARG";;
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

if [ -z "$panel" ]
then 
  # find target file for this batch from BatchInfo.txt with batch name.
  target_name=$(more /home/projects/cu_10184/projects/${dir_name}/Meta/BatchInfo.txt | awk -F '\t' -v batch="$batch" '( $1==batch ) {print $3}')
  if [ "${target_name}" = "" ]
  then 
    echo "Error: BatchInfo.txt does not have any information for this batch."
    exit 1
  fi
  target=/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/${target_name}.bed
elif [ "$panel" = "panel1" ]
then 
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38
  target=/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/${target_name}.bed
elif [ "$panel" = "panel2" ]
then 
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2
  target=/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/${target_name}.bed
else
  echo "Error: Please input panel version from command line as panel1 or panel2, or update BatchInfo.txt"
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
log_dir=${batch_dir}/log/ScanITD
rm -rf ${log_dir}
mkdir ${log_dir}

# error directory:
error_dir=${batch_dir}/error/ScanITD
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
# Lock for ITD:
Lock_ITD_dir=${Lock_dir}/ITD
Lock_ScanITD_dir=${Lock_ITD_dir}/ScanITD
rm -rf ${Lock_ScanITD_dir}
mkdir -p ${Lock_ScanITD_dir}

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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_ScanITD \
    -v batch=${batch},sample=${sample},BAM_dir=${BAM_dir},Lock_ScanITD_dir=${Lock_ScanITD_dir},target=${target},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Primary/ITD/ScanITD/ScanITD_job.sh
done

####################################################################################################################
####################################################################################################################
