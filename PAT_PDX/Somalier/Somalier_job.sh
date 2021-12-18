####################################################################################################################
####################################################################################################################
# Fingerprinting with Somalier - job.
# Author: Haiying Kong
# Last Modified: 26 October 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=10GB
#PBS -l walltime=20:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
sites="/home/projects/cu_10184/people/haikon/Software/somalier-0.2.13/sites.hg38.vcf.gz"

# Software tools:
Somalier="/home/projects/cu_10184/people/haikon/Software/somalier-0.2.13/somalier"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Somalier.
####################################################################################################################
${Somalier} extract -d ${dir_name}/Lock/Somalier/Extracted --sites ${sites} -f $hg ${pat_vcf}


####################################################################################################################
####################################################################################################################