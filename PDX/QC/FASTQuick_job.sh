####################################################################################################################
####################################################################################################################
# Quality control with QC - job.
# Author: Haiying Kong
# Last Modified: 22 October 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=40GB
#PBS -l walltime=200:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"

# Panel target:
target_nochr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Padded/${target_name}.bed
target_chr=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Padded/${target_name}.bed
target_nopad=/home/projects/cu_10184/projects/${dir_name}/PanelSeqData/Bait_Target/Target/Chr_Original/${target_name}.bed

# FASTQuick:
fastquick_reference_dir="/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37"
fastquick_resource_dir="/home/projects/cu_10184/people/haikon/Software/FASTQuick/resource"
target_fastquick=/home/projects/cu_10184/people/haikon/Reference/FASTQuick/hg37/${target_name}.bed

####################################################################################################################
# Software tools:
FASTQuick="/home/projects/cu_10184/people/haikon/Software/FASTQuick/bin/FASTQuick.sh"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Quality Control.
####################################################################################################################
# FASTQuick:
# Change to working directory.
cd ${temp_dir}

# Create a folder to save index files for the sample.
mkdir ${batch_dir}/QC/FASTQuick/Lock/${sample}

# Create a folder to save the results for the sample.
mkdir -p ${batch_dir}/QC/FASTQuick/Result/${sample}/tmp

# Run FASTQuick.
${FASTQuick}  \
  --steps All  \
  --fastq_1 ${fq_dir}/${sample}_R1.fq.gz  \
  --fastq_2 ${fq_dir}/${sample}_R2.fq.gz  \
  --output ${batch_dir}/QC/FASTQuick/Result/${sample}/${sample}  \
  --workingDir ${batch_dir}/QC/FASTQuick/Result/${sample}/tmp  \
  --reference ${fastquick_reference_dir}/hs37d5.fa  \
  --targetRegion ${target_fastquick}  \
  --dbSNP ${fastquick_reference_dir}/dbsnp132_20101103.vcf.gz  \
  --callableRegion ${fastquick_reference_dir}/20141020.strict_mask.whole_genome.bed  \
  --index ${batch_dir}/QC/FASTQuick/Lock/${sample}/index  \
  --candidateVCF ${fastquick_resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --SVDPrefix ${fastquick_resource_dir}/1000g.phase3.10k.b37.vcf.gz  \
  --nThread ${n_thread}

if [ $? -ne 0 ];  then echo "Failed at FASTQuick." >> ${batch_dir}/RedFlag/${sample}.txt;  fi

####################################################################################################################
####################################################################################################################
