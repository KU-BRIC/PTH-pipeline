####################################################################################################################
####################################################################################################################
# Run MuTect2 on tumor samples, and annotate variants.
# Author: Haiying Kong
# Last Modified: 8 August 2021
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

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}

# Change to working directory.
mkdir -p ${batch_dir}/temp
cd ${batch_dir}/temp

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/SNV_InDel/MuTect2_1
rm -rf ${log_dir}
log_dir=${log_dir}/Somatic
mkdir -p ${log_dir}

# error directory:
error_dir=${batch_dir}/error/SNV_InDel/MuTect2_1
rm -rf ${error_dir}
error_dir=${error_dir}/Somatic
mkdir -p ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
# Lock for variant calling:
rm -rf ${batch_dir}/Lock/SNV_InDel/MuTect2_1
mkdir -p ${batch_dir}/Lock/SNV_InDel/MuTect2_1
mkdir ${batch_dir}/Lock/SNV_InDel/MuTect2_1/vcf
mkdir ${batch_dir}/Lock/SNV_InDel/MuTect2_1/maf

####################################################################################################################
####################################################################################################################
# Define directories to save final results.
####################################################################################################################
rm -rf ${batch_dir}/Result/SNV_InDel/MuTect2_1
mkdir -p ${batch_dir}/Result/SNV_InDel/MuTect2_1

####################################################################################################################
####################################################################################################################
# Get BAM file names.
cd ${batch_dir}/Lock/BAM
bam_files=($(ls *.bam))

####################################################################################################################
# Change to working directory.
cd ${batch_dir}/temp

####################################################################################################################
# Get sample names.
samples=($(echo ${bam_files[@]%.bam} | tr ' ' '\n' | sort -u | tr '\n' ' '))

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
for sam in ${samples[@]}
do
  # If the sample is not NORMAL, submit job to run MuTect2.
  meta_file="/home/projects/cu_10184/projects/PTH/Meta/SampleInfo.txt"

  tail -n +2 ${meta_file} | awk -F"\t" -v b="${batch}" -v s="${sam}" '$1==b && $4==s && $5=="NORMAL"' >${batch_dir}/temp/${sam}.txt

  if ! [[ -s "${batch_dir}/temp/${sam}.txt" ]]
  then
    qsub -o ${log_dir}/${sam}.log -e ${error_dir}/${sam}.error -N ${batch}_${sam}_MuTect2_1  \
      -v target_name=${target_name},dir_name=${dir_name},batch_dir=${batch_dir},sam=${sam}  \
      /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/Somatic/MuTect2_job.sh
  fi

  rm ${batch_dir}/temp/${sam}.txt
done

####################################################################################################################
####################################################################################################################
