####################################################################################################################
####################################################################################################################
# Pre-process, call variants, annotate and filter variants.
# Author: Haiying Kong
# Last Modified: 12 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

source /home/projects/cu_10184/projects/PTH/Software/envsetup

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":d:b:t:c:h:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    b) batch="$OPTARG";;
    t) n_thread="$OPTARG";;
    c) classify_tool="$OPTARG";;
    h) thresh_n_human_read="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&5;;
  esac
done
n_thread=${n_thread#0}
thresh_n_human_read=${thresh_n_human_read#0}

####################################################################################################################
####################################################################################################################
# Check if arguments are input from command.
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

if [ -z "${n_thread}" ]
then
  n_thread=8
fi

####################################################################################################################
####################################################################################################################
# Define directory for the batch.
fq_dir=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/${batch}/fastq/${classify_tool}
batch_dir=/home/projects/cu_10184/projects/${dir_name}/BatchWork/${batch}

# Create a new directory for the batch.
mkdir -p ${batch_dir}

# Change to working directory.
temp_dir=${batch_dir}/temp
mkdir -p ${temp_dir}
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save log files and error files.
# log directory:
mkdir -p ${batch_dir}/log
log_dir=${batch_dir}/log/Ensemble
rm -rf ${log_dir}
mkdir ${log_dir}

# error directory:
mkdir -p ${batch_dir}/error
error_dir=${batch_dir}/error/Ensemble
rm -rf ${error_dir}
mkdir ${error_dir}

####################################################################################################################
####################################################################################################################
# Define directories to save intermediate results.
####################################################################################################################
Lock_dir=${batch_dir}/Lock
rm -rf ${Lock_dir}
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
rm -rf ${Result_dir}
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
# find target file for this batch.
target_name=$(more /home/projects/cu_10184/projects/${dir_name}/Meta/BatchInfo.txt | awk -F '\t' -v batch="$batch" '( $1==batch ) {print $3}')
target_nochr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
target_chr=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed
target_nopad=/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Original/${target_name}.bed

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
# samples=$(awk -F"\t" -v var=${pct_thresh} 'NR>1 && $8>var' /home/projects/cu_10184/projects/PTH/PanelSeqData/PDX_001/meta/Summary_ReadCounts.txt | cut -d$'\t' -f 1)

####################################################################################################################
# Run pipeline on all samples in this batch.
####################################################################################################################
for sample in ${samples[@]}
do
  # Check if the sample has sufficient reads from human genome.
  read_count_file=/home/projects/cu_10184/projects/${dir_name}/PDXseqData/${batch}/${classify_tool}/ReadCount/${sample}.txt
  n_human_read=$(awk -F"\t" 'NR>1 {print $2}' ${read_count_file})
  if [ ${n_human_read} -lt ${thresh_n_human_read} ]
    then
      echo "$batch $sample: Too few reads from human genome."
    else
      qsub -o ${log_dir}/${sample}.log -e ${error_dir}/${sample}.error -N ${batch}_${sample}_Ensemble \
        -v n_thread=${n_thread},target_chr=${target_chr},target_nochr=${target_nochr},target_nopad=${target_nopad},batch=${batch},sample=${sample},fq_dir=${fq_dir},BAM_dir=${BAM_dir},BAM_lock_dir=${BAM_lock_dir},Lock_SNV_InDel_dir=${Lock_SNV_InDel_dir},Lock_VarDict_dir=${Lock_VarDict_dir},Lock_SNVer_dir=${Lock_SNVer_dir},Lock_LoFreq_dir=${Lock_LoFreq_dir},Result_SNV_InDel_dir=${Result_SNV_InDel_dir},Lock_CNACS_dir=${Lock_CNACS_dir},Lock_DOC_dir=${Lock_DOC_dir},Lock_SoftClipping_dir=${Lock_SoftClipping_dir},temp_dir=${temp_dir} \
        /home/projects/cu_10184/projects/PTH/Code/PDX/Ensemble/Ensemble_job.sh
  fi
done

####################################################################################################################
####################################################################################################################
