####################################################################################################################
####################################################################################################################
# Run CNACS for all samples in the batch - job.
# Author: Haiying Kong
# Last Modified: 25 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=70GB
#PBS -l walltime=200:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Call CNVs with CNACS.
####################################################################################################################
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
