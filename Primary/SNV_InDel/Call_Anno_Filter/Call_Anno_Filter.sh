####################################################################################################################
####################################################################################################################
# Call variants, annotate and filter variants.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 8 August 2021
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
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/SNV_InDel
rm -rf ${log_dir}
mkdir ${log_dir}

# error directory:
error_dir=${batch_dir}/error/SNV_InDel
rm -rf ${error_dir}
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
# Lock for variant calling:
rm -rf ${batch_dir}/Lock/SNV_InDel
mkdir ${batch_dir}/Lock/SNV_InDel
mkdir ${batch_dir}/Lock/SNV_InDel/VarDict
mkdir ${batch_dir}/Lock/SNV_InDel/VarDict/vcf_0
mkdir ${batch_dir}/Lock/SNV_InDel/VarDict/vcf
mkdir ${batch_dir}/Lock/SNV_InDel/VarDict/maf
mkdir ${batch_dir}/Lock/SNV_InDel/SNVer
mkdir ${batch_dir}/Lock/SNV_InDel/SNVer/vcf
mkdir ${batch_dir}/Lock/SNV_InDel/SNVer/maf
mkdir ${batch_dir}/Lock/SNV_InDel/LoFreq
mkdir ${batch_dir}/Lock/SNV_InDel/LoFreq/vcf
mkdir ${batch_dir}/Lock/SNV_InDel/LoFreq/maf

####################################################################################################################
####################################################################################################################
# Define directories to save final results.
####################################################################################################################
mkdir -p ${batch_dir}/Result

####################################################################################################################
# Result for variant calling:
rm -rf ${batch_dir}/Result/SNV_InDel
mkdir ${batch_dir}/Result/SNV_InDel
mkdir ${batch_dir}/Result/SNV_InDel/AllVariants
mkdir ${batch_dir}/Result/SNV_InDel/AllVariants/Callers_Long
mkdir ${batch_dir}/Result/SNV_InDel/AllVariants/Callers_Wide
mkdir ${batch_dir}/Result/SNV_InDel/Filtered
mkdir ${batch_dir}/Result/SNV_InDel/Filtered/Long
mkdir ${batch_dir}/Result/SNV_InDel/Filtered/Medium
mkdir ${batch_dir}/Result/SNV_InDel/Filtered/Short

####################################################################################################################
####################################################################################################################
# Get BAM file names.
cd ${batch_dir}/Lock/BAM
bam_files=($(ls *.bam))

####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Get sample names.
samples=($(echo ${bam_files[@]%.bam} | tr ' ' '\n' | sort -u | tr '\n' ' '))

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
for sample in ${samples[@]}
do
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_SNV_InDel \
    -v n_thread=${n_thread},target_name=${target_name},dir_name=${dir_name},batch=${batch},sample=${sample},batch_dir=${batch_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Call_Anno_Filter/Call_Anno_Filter_job.sh
done

####################################################################################################################
####################################################################################################################
