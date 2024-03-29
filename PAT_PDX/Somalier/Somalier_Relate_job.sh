####################################################################################################################
####################################################################################################################
# Fingerprinting with Somalier to relate - job.
# Author: Haiying Kong
# Last Modified: 29 October 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=10GB
#PBS -l walltime=20:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
# Software tools:
Somalier="/home/projects/cu_10184/people/haikon/Software/somalier-0.2.13/somalier"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to result directory.
cd ${res_dir}

####################################################################################################################
# Run Somalier.
shopt -s globstar
${Somalier} relate ${lock_dir}/**/*.somalier

####################################################################################################################
####################################################################################################################
