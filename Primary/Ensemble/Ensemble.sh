####################################################################################################################
####################################################################################################################
# Pre-process, call variants, annotate and filter variants.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 16 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

source /home/projects/cu_10184/projects/PTH/Software/envsetup

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
    \?) echo "Invalid option -$OPTARG" >&3;;
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

# Create a new directory for the batch.
rm -rf ${batch_dir}
mkdir -p ${batch_dir}

# Change to working directory.
temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
mkdir ${batch_dir}/log
log_dir=${batch_dir}/log/Ensemble
mkdir ${log_dir}

# error directory:
mkdir ${batch_dir}/error
error_dir=${batch_dir}/error/Ensemble
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
Lock_dir=${batch_dir}/Lock
mkdir ${Lock_dir}

####################################################################################################################
# Lock for BAM:
BAM_dir=${Lock_dir}/BAM
mkdir ${BAM_dir}

BAM_lock_dir=${BAM_dir}/lock
mkdir ${BAM_lock_dir}

####################################################################################################################
# Lock for variant calling:
Lock_SNV_InDel_dir=${Lock_dir}/SNV_InDel
mkdir ${Lock_SNV_InDel_dir}

Lock_VarDict_dir=${Lock_SNV_InDel_dir}/VarDict
mkdir ${Lock_VarDict_dir}
mkdir ${Lock_VarDict_dir}/vcf
mkdir ${Lock_VarDict_dir}/maf

Lock_SNVer_dir=${Lock_SNV_InDel_dir}/SNVer
mkdir ${Lock_SNVer_dir}
mkdir ${Lock_SNVer_dir}/vcf
mkdir ${Lock_SNVer_dir}/maf

Lock_LoFreq_dir=${Lock_SNV_InDel_dir}/LoFreq
mkdir ${Lock_LoFreq_dir}
mkdir ${Lock_LoFreq_dir}/vcf
mkdir ${Lock_LoFreq_dir}/maf

####################################################################################################################
# Lock for CNV:
Lock_CNV_dir=${Lock_dir}/CNV
mkdir ${Lock_CNV_dir}
Lock_CNACS_dir=${Lock_CNV_dir}/CNACS
mkdir ${Lock_CNACS_dir}

####################################################################################################################
# Lock for depth of coverage:
Lock_DOC_dir=${Lock_dir}/DepthOfCoverage
mkdir ${Lock_DOC_dir}
mkdir ${Lock_DOC_dir}/FreqTable
mkdir ${Lock_DOC_dir}/DensityPlot

####################################################################################################################
# Lock for ITD:
Lock_ITD_dir=${Lock_dir}/ITD
Lock_SoftClipping_dir=${Lock_ITD_dir}/SoftClipping

####################################################################################################################
####################################################################################################################
# Define directories to save final results.
####################################################################################################################
Result_dir=${batch_dir}/Result
mkdir ${Result_dir}

####################################################################################################################
# Result for variant calling:
Result_SNV_InDel_dir=${Result_dir}/SNV_InDel
mkdir ${Result_SNV_InDel_dir}
mkdir ${Result_SNV_InDel_dir}/AllVariants
mkdir ${Result_SNV_InDel_dir}/AllVariants/Callers_Long
mkdir ${Result_SNV_InDel_dir}/AllVariants/Callers_Wide
mkdir ${Result_SNV_InDel_dir}/Filtered
mkdir ${Result_SNV_InDel_dir}/Filtered/Long
mkdir ${Result_SNV_InDel_dir}/Filtered/Medium
mkdir ${Result_SNV_InDel_dir}/Filtered/Short

Result_CNV_dir=${Result_dir}/CNV
mkdir ${Result_CNV_dir}
Result_CNACS_dir=${Result_CNV_dir}/CNACS
mkdir ${Result_CNACS_dir}

Result_DOC_dir=${Result_dir}/DOC
mkdir ${Result_DOC_dir}

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
  qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_Ensemble \
    -v n_thread=${n_thread},target_chr=${target_chr},target_nochr=${target_nochr},target_nopad=${target_nopad},batch=${batch},sample=${sample},fq_dir=${fq_dir},BAM_dir=${BAM_dir},BAM_lock_dir=${BAM_lock_dir},Lock_SNV_InDel_dir=${Lock_SNV_InDel_dir},Lock_VarDict_dir=${Lock_VarDict_dir},Lock_SNVer_dir=${Lock_SNVer_dir},Lock_LoFreq_dir=${Lock_LoFreq_dir},Result_SNV_InDel_dir=${Result_SNV_InDel_dir},Lock_CNACS_dir=${Lock_CNACS_dir},Lock_DOC_dir=${Lock_DOC_dir},Lock_SoftClipping_dir=${Lock_SoftClipping_dir},temp_dir=${temp_dir} \
    /home/projects/cu_10184/projects/PTH/Code/Primary/Ensemble/Ensemble_job.sh
done

####################################################################################################################
####################################################################################################################
