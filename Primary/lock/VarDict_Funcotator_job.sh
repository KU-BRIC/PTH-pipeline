####################################################################################################################
####################################################################################################################
# Run VarDict and annotate with Funcotator - job.
# Author: Haiying Kong
# Last Modified: 1 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=70GB
#PBS -l walltime=200:00:00

source /home/people/haikon/.bashrc

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10145/people/haikon/Reference/GATK/hg38/Homo_sapiens_assembly38.fasta"
Funcotator_DB="/home/projects/cu_10145/people/haikon/Reference/Funcotator/funcotator_dataSources.v1.7.20200521s"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Call SNV_InDel variants and annotate.
####################################################################################################################
# VarDict:
################################################
# Set parameters.
allel_freq=0.01
# map_qual=25
phred=22.5
var_pos=5
Qratio=1.5
indel_size=50
min_match=25
min_SV=500
ins_size=250
ins_SD=100
nmfreq=0.1
mfreq=0.25

# Run VarDict:
vardict -G $hg -N $sample -b ${BAM_dir}/${sample}.bam $target_nochr -c 1 -S 2 -E 3 \
  -f ${allel_freq} -q ${phred} -P ${var_pos} -o ${Qratio} -I ${indel_size}  \
  -M ${min_match} -L ${min_SV} -w ${ins_size} -W ${ins_SD} --nmfreq ${nmfreq} --mfreq ${mfreq}  \
  | teststrandbias.R \
  | var2vcf_valid.pl -N ${sample} -E -f ${allel_freq} > ${Lock_VarDict_dir}/vcf/${sample}.vcf
# -O ${map_qual}

# Annotation with Funcotator:
gatk Funcotator \
  -R $hg -V ${Lock_VarDict_dir}/vcf/$sample.vcf -O ${Lock_VarDict_dir}/maf/${sample}.maf \
  --output-file-format MAF --ref-version hg38 --data-sources-path ${Funcotator_DB} \
  --disable-sequence-dictionary-validation

# Fill missing DP and AF in Funcotator maf output.
Rscript /home/projects/cu_10184/projects/PTH/Code/Ensemble/Fill_DP_AF_FuncotatorMAF_nonSNVer_job.R ${Lock_VarDict_dir} $sample 'VarDict'

####################################################################################################################
####################################################################################################################
