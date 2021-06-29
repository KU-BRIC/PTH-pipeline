####################################################################################################################
####################################################################################################################
# Perform fingerprinting with verifyBamID.
# Author: Balthasar Schlotmann
# Last Modified: 10 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:t:r:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    t) n_thread="$OPTARG";;
    r) run_name="$OPTARG";;
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
log_dir=${qc_dir}/log/VerifyBamID/${run_name}
rm -rf ${log_dir}
mkdir -p ${log_dir}/verifyBamID

# error directory:
error_dir=${qc_dir}/error/VerifyBamID/${run_name}
rm -rf ${error_dir}
mkdir -p ${error_dir}

# intermediate results:
lock_dir=${qc_dir}/Lock/VerifyBamID/${run_name}
rm -rf ${lock_dir}
mkdir -p ${lock_dir}

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
# Read in the sample pair list line by line and submit job for fingerprinting.
sed 1d /home/projects/cu_10184/projects/${dir_name}/QC/Meta/${run_name}SampleList.txt | while read f
do

  batch1=$(echo $f | cut -f1 -d' ')
  sample1=$(echo $f | cut -f2 -d' ')
  batch2=$(echo $f | cut -f3 -d' ')
  sample2=$(echo $f | cut -f4 -d' ')

  vcf=${qc_dir}/Lock/VerifyBamID/Genotype/${batch1}/${sample1}.vcf
  bam=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch2}/Lock/BAM/${sample2}.bam

  echo "  qsub -o ${log_dir}/${sample1}_vs_${sample2}.log -e ${error_dir}/${sample1}_vs_${sample2}.error -N ${batch1}_${sample1}_vs_${sample2}_VerifyBamID  \
    -v n_thread=${n_thread},sample1=${sample1},sample2=${sample2},vcf=$vcf,bam=${bam},lock_dir=${lock_dir},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/QC/VerifyBamID/VerifyBamID_job.sh"


  qsub -o ${log_dir}/${sample1}_vs_${sample2}.log -e ${error_dir}/${sample1}_vs_${sample2}.error -N ${batch1}_${sample1}_vs_${sample2}_VerifyBamID  \
    -v n_thread=${n_thread},sample1=${sample1},sample2=${sample2},vcf=$vcf,bam=${bam},lock_dir=${lock_dir},log_dir=${log_dir},temp_dir=${temp_dir}  \
    /home/projects/cu_10184/projects/PTH/Code/QC/VerifyBamID/VerifyBamID_job.sh

done

####################################################################################################################
####################################################################################################################
