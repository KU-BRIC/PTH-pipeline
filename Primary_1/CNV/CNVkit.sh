####################################################################################################################
####################################################################################################################
# Run CNVkit for all samples in the batch.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 5 September 2021
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
temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/CNVkit
rm -rf ${log_dir}
mkdir ${log_dir}

# error directory:
error_dir=${batch_dir}/error/CNVkit
rm -rf ${error_dir}
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
# Lock for CNV:
rm -rf ${batch_dir}/Lock/CNV/CNVkit
mkdir -p ${batch_dir}/Lock/CNV/CNVkit

# Result for CNV:
rm -rf ${batch_dir}/Result/CNV/CNVkit
mkdir -p ${batch_dir}/Result/CNV/CNVkit/Annotate
mkdir -p ${batch_dir}/Result/CNV/CNVkit/XCNV
mkdir -p ${batch_dir}/Result/CNV/CNVkit/AnnotSV

####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
meta_file="/home/projects/cu_10184/projects/PTH/Meta/SampleInfo.txt"

for sam in ${batch_dir}/Lock/BAM/*.bam
do
  sam=$(basename $sam)
  sam=${sam/.bam/}

  # If the sample is not NORMAL, submit job to run MuTect2.
  tail -n +2 ${meta_file} | awk -F"\t" -v b="${batch}" -v s="${sam}" '$1==b && $4==s && $5=="NORMAL"' >${batch_dir}/temp/${sam}.txt

  if ! [[ -s "${batch_dir}/temp/${sam}.txt" ]]
  then
    qsub -o ${log_dir}/${sam}.log -e ${error_dir}/${sam}.error -N ${batch}_${sam}_CNVkit \
      -v sam=${sam},batch_dir=${batch_dir},temp_dir=${temp_dir} \
      /home/projects/cu_10184/projects/PTH/Code/Primary_1/CNV/CNVkit_job.sh
  fi
  rm ${batch_dir}/temp/${sam}.txt
done

####################################################################################################################
####################################################################################################################
