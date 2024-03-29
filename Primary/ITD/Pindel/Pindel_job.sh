####################################################################################################################
####################################################################################################################
# Identify ITD with Pindel - job.
# Author: Haiying Kong
# Last Modified: 18 April 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=10:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"

# Set parameters.
n_thread=8

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run Pindel:
echo -e ${BAM_dir}/${sample}.bam"\t"250"\t"${sample} > ${Lock_Pindel_dir}/${sample}_config.txt
pindel -f $hg -t $target -c chr13:28033736-28034557  \
       -i ${Lock_Pindel_dir}/${sample}_config.txt -o ${Lock_Pindel_dir}/$sample  \
       -T $n_thread -x 3 -u 0.04 
rm ${Lock_Pindel_dir}/${sample}_config.txt

pindel2vcf -p ${Lock_Pindel_dir}/${sample}_TD -r $hg -R HG38 -d 20201224 -v ${Lock_Pindel_dir}/${sample}.vcf

conda deactivate

####################################################################################################################
####################################################################################################################
