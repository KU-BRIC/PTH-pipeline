####################################################################################################################
####################################################################################################################
# Run LoFreq.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 31 May 2021
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
  target_nochr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
  target_chr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed
  target_nopad=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Original/${target_name}.bed
elif [ "$panel" = "panel1" ]
then
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38
  target_nochr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
  target_chr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed
  target_nopad=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Original/${target_name}.bed
elif [ "$panel" = "panel2" ]
then
  # Find target file for panel version 1.
  target_name=all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2
  target_nochr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
  target_chr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed
  target_nopad=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Original/${target_name}.bed
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
fq_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

# Change to working directory.
temp_dir=${batch_dir}/temp
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
log_dir=${batch_dir}/log/LoFreq
rm -rf ${log_dir}
mkdir ${log_dir}

# error directory:
error_dir=${batch_dir}/error/LoFreq
rm -rf ${error_dir}
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
Lock_dir=${batch_dir}/Lock
BAM_dir=${Lock_dir}/BAM

# Lock for variant calling:
Lock_SNV_InDel_dir=${Lock_dir}/SNV_InDel
Lock_LoFreq_dir=${Lock_SNV_InDel_dir}/LoFreq
rm -rf ${Lock_LoFreq_dir}
mkdir ${Lock_LoFreq_dir}
mkdir ${Lock_LoFreq_dir}/vcf
mkdir ${Lock_LoFreq_dir}/maf

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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_LoFreq \
    -v n_thread=${n_thread},target_chr=${target_chr},target_nochr=${target_nochr},target_nopad=${target_nopad},batch=${batch},sample=${sample},BAM_dir=${BAM_dir},Lock_LoFreq_dir=${Lock_LoFreq_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/LoFreq_job.sh
done

####################################################################################################################
####################################################################################################################
