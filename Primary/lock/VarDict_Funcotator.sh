####################################################################################################################
####################################################################################################################
# Run VarDict and annotate with Funcotator.
# Author: Haiying Kong
# Last Modified: 1 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

source /home/people/haikon/.bashrc

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":b:t:" opt
do
  case $opt in
    b) batch="$OPTARG";;
    t) n_thread="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done
n_thread=${n_thread#0}

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
fq_dir=/home/projects/cu_10184/projects/PTH/PanelSeqData/${batch}/fastq
batch_dir=/home/projects/cu_10184/projects/PTH/BatchWork/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
if [ ! -d "${temp_dir}" ]
  then mkdir ${temp_dir}
fi
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/VarDict
if [ -d "${log_dir}" ]
  then rm -r ${log_dir}
fi
mkdir ${log_dir}

# error directory:
error_dir=${batch_dir}/error/VarDict
if [ -d "${error_dir}" ]
  then rm -r ${error_dir}
fi
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
Lock_dir=${batch_dir}/Lock

BAM_dir=${Lock_dir}/BAM

####################################################################################################################
# Lock for variant calling:
Lock_SNV_InDel_dir=${Lock_dir}/SNV_InDel
Lock_VarDict_dir=${Lock_SNV_InDel_dir}/VarDict
if [ -d "${Lock_VarDict_dir}" ]
  then rm -r ${Lock_VarDict_dir}
fi
mkdir ${Lock_VarDict_dir}
mkdir ${Lock_VarDict_dir}/vcf
mkdir ${Lock_VarDict_dir}/maf

####################################################################################################################
####################################################################################################################
# find target file for this batch.
target_name=$(more /home/projects/cu_10184/projects/PTH/Meta/BatchInfo.txt | awk -F '\t' -v batch="$batch" '( $1==batch ) {print $3}')
target_nochr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
target_chr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed

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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_VarDict \
    -v n_thread=${n_thread},target_chr=${target_chr},target_nochr=${target_nochr},sample=${sample},BAM_dir=${BAM_dir},Lock_SNV_InDel_dir=${Lock_SNV_InDel_dir},Lock_VarDict_dir=${Lock_VarDict_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Ensemble/lock/VarDict_Funcotator_job.sh
done

####################################################################################################################
####################################################################################################################
