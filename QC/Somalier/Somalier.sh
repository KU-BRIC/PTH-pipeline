####################################################################################################################
####################################################################################################################
# Perform sample quanlity control with FASTQuick.
# Author: Balthasar Schlotmann
# Last Modified: 24 June 2021
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

####################################################################################################################
####################################################################################################################
# Define directories.
qc_dir=/home/projects/cu_10184/projects/${dir_name}/QC

# Change to working directory.
temp_dir=/home/projects/cu_10184/projects/${dir_name}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

# log directory:
log_dir=${qc_dir}/log/Somalier
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${qc_dir}/error/Somalier
rm -rf ${error_dir}
mkdir -p ${error_dir}

# fastq directory:
fq_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq

# intermediate file directory:
lock_dir=${qc_dir}/Lock/Somalier/${batch}
rm -rf ${lock_dir}
mkdir -p ${lock_dir}

####################################################################################################################
# Lock for BAM:
BAM_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}/Lock/BAM/

####################################################################################################################
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

shopt -s nullglob

for sample in ${BAM_dir}/*.bam
do
  bam=$sample
  sample=$(basename $sample)
  sample=${sample/.bam/}
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_Somalier  \
    -v bam=${bam},batch=${batch},sample=${sample},lock_dir=${lock_dir},res_dir=${res_dir},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/QC/Somalier/Somalier_job.sh
done

####################################################################################################################
####################################################################################################################
