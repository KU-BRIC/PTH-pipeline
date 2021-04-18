####################################################################################################################
####################################################################################################################
# Identify ITD with ScanITD - job.
# Author: Haiying Kong
# Last Modified: 18 April 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=10:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate
conda activate ScanITD

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"

# Software tools:
ScanITD="/home/projects/cu_10184/people/haikon/Software/ScanITD/ScanITD.py"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run ScanITD:
$ScanITD -r $hg -t $target -i ${BAM_dir}/${sample}.bam -o ${Lock_ScanITD_dir}/${sample} -m 20 -c 1 -d 100 -f 0.0001 -l 5 -n 3
conda deactivate

#  -m MAPQ, --mapq MAPQ  minimal MAPQ in BAM for calling ITD (default: 15)
#  -c AO, --ao AO        minimal observation count for ITD (default: 4)
#  -d DP, --depth DP     minimal depth to call ITD (default: 10)
#  -f VAF, --vaf VAF     minimal variant allele frequency (default: 0.1)
#  -l ITD_LEN, --len ITD_LEN
#                        minimal ITD length to report (default: 10)
#  -n MISMATCH           maximum allowed mismatch bases of pairwise local alignment (default: 3)

####################################################################################################################
####################################################################################################################
