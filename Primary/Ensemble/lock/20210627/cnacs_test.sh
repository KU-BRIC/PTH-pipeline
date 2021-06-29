####################################################################################################################
####################################################################################################################
# Pre-process, call variants, annotate and filter variants.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 23 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch name and number of thread from argument passing.
dir_name=PTH_test
batch=ToyBatch
panel=panel2
n_thread=8

####################################################################################################################
####################################################################################################################
# Check if argument inputs from command are correct.
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
# Define directories to save intermediate results.
####################################################################################################################
Lock_dir=${batch_dir}/Lock

####################################################################################################################
# Lock for BAM:
BAM_dir=${Lock_dir}/BAM

####################################################################################################################
# Lock for CNV:
Lock_CNV_dir=${Lock_dir}/CNV
Lock_CNACS_dir=${Lock_CNV_dir}/CNACS

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
sample=${samples[0]}

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
Funcotator_DB="/home/projects/cu_10184/people/haikon/Reference/Funcotator/funcotator_dataSources.v1.7.20200521s"

conda activate CNACS

export JAVAPATH=/home/projects/cu_10184/people/haikon/Software/anaconda3/bin
export PICARD_PATH=/home/projects/cu_10184/people/haikon/Software/picard.jar
export SAMTOOLS_PATH=/home/projects/cu_10184/people/haikon/Software/samtools-1.12/bin
export PERL_PATH=/usr/bin/perl
export BEDTOOLS_PATH=/home/projects/cu_10184/people/haikon/Software/anaconda3/bin
export R_PATH=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/R
export R_LIBS_PATH=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/library
export R_LIBS=/home/projects/cu_10184/people/haikon/Software/R-4.0.4/library

rm -rf ${Lock_CNACS_dir}/${sample}
mkdir -p ${Lock_CNACS_dir}/${sample}

toil_cnacs run \
    ${Lock_CNACS_dir}/${sample}/jobstore/ \
    --stats \
    --db_dir /home/projects/cu_10184/projects/PTH/Reference/CNACS/db \
    --outdir ${Lock_CNACS_dir} \
    --pool_dir /home/projects/cu_10184/projects/PTH/Reference/CNACS/PoN \
    --fasta $hg \
    --samp ${BAM_dir}/${sample}.bam

conda deactivate


####################################################################################################################
####################################################################################################################
