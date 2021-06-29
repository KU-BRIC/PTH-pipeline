####################################################################################################################
####################################################################################################################
# Run CNVkit for all samples in the batch - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 22 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=40GB
#PBS -l walltime=20:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
output_reference="/home/projects/cu_10184/projects/PTH/Reference/CNVkit/OutputReference.cnn"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
####################################################################################################################
# Call CNVs with CNVkit.
####################################################################################################################
conda activate CNVkit

cnvkit.py batch ${BAM_dir}/${sample}.bam \
    --reference ${output_reference} \
    --output-dir ${Lock_CNVkit_dir}

echo "cnvkit.py scatter ${Lock_CNVkit_dir}/${sample}.cnr -s ${Lock_CNVkit_dir}/${sample}.cns -o ${Lock_CNVkit_dir}/${sample}.scatter.pdf"

cnvkit.py scatter ${Lock_CNVkit_dir}/${sample}.cnr -s ${Lock_CNVkit_dir}/${sample}.cns -o ${Lock_CNVkit_dir}/${sample}.scatter.pdf
cnvkit.py diagram ${Lock_CNVkit_dir}/${sample}.cnr -s ${Lock_CNVkit_dir}/${sample}.cns -o ${Lock_CNVkit_dir}/${sample}.diagram.pdf

conda deactivate

####################################################################################################################
####################################################################################################################
