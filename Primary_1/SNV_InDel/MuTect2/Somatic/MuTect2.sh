####################################################################################################################
####################################################################################################################
# Run MuTect2 on tumor samples, and annotate variants.
# Author: Haiying Kong
# Last Modified: 21 July 2021
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
# Define directory for the batch.
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}

# Change to working directory.
mkdir -p ${batch_dir}/temp
cd ${batch_dir}/temp

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/SNV_InDel/MuTect2
rm -rf ${log_dir}
log_dir=${log_dir}/Somatic
mkdir -p ${log_dir}

# error directory:
error_dir=${batch_dir}/error/SNV_InDel/MuTect2
rm -rf ${error_dir}
error_dir=${error_dir}/Somatic
mkdir -p ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
# Lock for variant calling:
rm -rf ${batch_dir}/Lock/SNV_InDel/MuTect2
mkdir -p ${batch_dir}/Lock/SNV_InDel/MuTect2
mkdir ${batch_dir}/Lock/SNV_InDel/MuTect2/vcf
mkdir ${batch_dir}/Lock/SNV_InDel/MuTect2/maf

####################################################################################################################
####################################################################################################################
# Define directories to save final results.
####################################################################################################################
rm -rf ${batch_dir}/Result/SNV_InDel/MuTect2
mkdir -p ${batch_dir}/Result/SNV_InDel/MuTect2

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

  meta_file=/home/projects/cu_10184/projects/${dir_name}/Meta/SampleInfo.txt
  tail -n +2 ${meta_file} | awk -F"\t" -v b="$batch" -v s="${sam}" '$1==b && $4==s && $5=="NORMAL"' >${batch_dir}/temp/${sam}.txt

  if ! [[ -s "${batch_dir}/temp/${sam}.txt" ]]
  then
    qsub -o ${log_dir}/${sam}.log -e ${error_dir}/${sam}.error -N ${batch}_${sam}_MuTect2  \
      -v sam=${sam},batch_dir=${batch_dir}  \
      /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2/Somatic/MuTect2_job.sh
  fi

  rm ${batch_dir}/temp/${sam}.txt
done

####################################################################################################################
####################################################################################################################
