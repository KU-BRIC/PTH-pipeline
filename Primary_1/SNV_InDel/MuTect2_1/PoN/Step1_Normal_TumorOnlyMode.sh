####################################################################################################################
####################################################################################################################
# Run MuTect2 on panel2 normals in tumor only mode.
# Author: Haiying Kong
# Last Modified: 22 July 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done

####################################################################################################################
####################################################################################################################
# Check if argument inputs from command are correct.
if [[ -z "${dir_name}" ]]
then
  echo "By default, directory name is set to PTH."
  dir_name=PTH
fi

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
all_dir=/home/projects/cu_10184/projects/${dir_name}/AllBatches_1
mkdir -p ${all_dir}

# Change to working directory.
temp_dir=${all_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${all_dir}/log/SNV_InDel/MuTect2_1/PoN_Step1
rm -rf ${log_dir}
mkdir -p ${log_dir}

# error directory:
error_dir=${all_dir}/error/SNV_InDel/MuTect2_1/PoN_Step1
rm -rf ${error_dir}
mkdir -p ${error_dir}

####################################################################################################################
# Create directory to save results.
rm -rf ${all_dir}/Lock/SNV_InDel/MuTect2_1/PoN/Normal
mkdir -p ${all_dir}/Lock/SNV_InDel/MuTect2_1/PoN/Normal

####################################################################################################################
####################################################################################################################
# Read in meta file line by line and if it is normal sample, submit job to run MuTect2.
meta_file="/home/projects/cu_10184/projects/PTH/Meta/SampleInfo_Normal_Panel2.txt"

tail -n +2 ${meta_file} | while IFS= read -r one_row
do

  # Get batch and sample name.
  tmp=($(echo "${one_row}" | tr '\t' '\n'))
  batch=${tmp[0]}
  sam=${tmp[3]}
    
  # Input and output files.
  bam=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Lock/BAM/${sam}.bam
  vcf=${all_dir}/Lock/SNV_InDel/MuTect2_1/PoN/Normal/${sam}.vcf.gz

  # Submit job.
  qsub -o ${log_dir}/${sam}.log -e ${error_dir}/${sam}.error -N ${batch}_${sam}_MuTect2_1_PoN_Step1  \
    -v bam=${bam},vcf=${vcf},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/PoN/Step1_Normal_TumorOnlyMode_job.sh

done

####################################################################################################################
####################################################################################################################
