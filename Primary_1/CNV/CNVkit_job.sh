####################################################################################################################
####################################################################################################################
# Run CNVkit for all samples in the batch - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 5 September 2021
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
refFlat="/home/projects/cu_10184/people/haikon/Reference/UCSC/hg38/refFlat.txt"

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

cnvkit.py batch ${batch_dir}/Lock/BAM/${sam}.bam \
    --reference ${output_reference}  \
    --output-dir ${batch_dir}/Lock/CNV/CNVkit

echo "cnvkit.py scatter ${batch_dir}/Lock/CNV/CNVkit/${sam}.cnr -s ${batch_dir}/Lock/CNV/CNVkit/${samp}.cns -o ${batch_dir}/Lock/CNV/CNVkit/${sam}.scatter.pdf"

cnvkit.py scatter ${batch_dir}/Lock/CNV/CNVkit/${sam}.cnr -s ${batch_dir}/Lock/CNV/CNVkit/${sam}.cns -o ${batch_dir}/Lock/CNV/CNVkit/${sam}.scatter.pdf
cnvkit.py diagram ${batch_dir}/Lock/CNV/CNVkit/${sam}.cnr -s ${batch_dir}/Lock/CNV/CNVkit/${sam}.cns -o ${batch_dir}/Lock/CNV/CNVkit/${sam}.diagram.pdf
cnv_annotate.py ${refFlat} ${batch_dir}/Lock/CNV/CNVkit/${sam}.cns -o ${batch_dir}/Result/CNV/CNVkit/Annotate/${sam}.cns

conda deactivate

####################################################################################################################


####################################################################################################################
####################################################################################################################
