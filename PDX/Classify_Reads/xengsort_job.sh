####################################################################################################################
####################################################################################################################
# Classify xenofraft sample reads to human and mouse - job.
# Author: Haiying Kong
# Last Modified: 25 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=50GB
#PBS -l walltime=200:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
mouse_human_index="/home/projects/cu_10184/people/haikon/Reference/UCSC/xengsort/Index_hg38_mm39.h5"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# Run xengsort to classify the reads.
conda activate xengsort

gunzip -c ${in_dir}/${sample}_R1.fq.gz >${out_dir}/${sample}_R1.fq
gunzip -c ${in_dir}/${sample}_R2.fq.gz >${out_dir}/${sample}_R2.fq

xengsort classify -T ${n_thread} --index ${mouse_human_index}  \
  --fastq ${out_dir}/${sample}_R1.fq --pairs ${out_dir}/${sample}_R2.fq  \
  --prefix ${out_dir}/${sample}

conda deactivate

# Save read counts for reads in all categories.
echo "# All Human Mouse Both Neither Ambiguous" > ${out_dir}/ReadCount/${sample}.txt
echo $(expr $(wc -l < ${out_dir}/${sample}_R1.fq) / 4)$'\t'$(expr $(wc -l < ${out_dir}/${sample}-graft.1.fq) / 4)$'\t'$(expr $(wc -l < ${out_dir}/${sample}-host.1.fq) / 4)$'\t'$(expr $(wc -l < ${out_dir}/${sample}-both.1.fq) / 4)$'\t'$(expr $(wc -l < ${out_dir}/${sample}-neither.1.fq) / 4)$'\t'$(expr $(wc -l < ${out_dir}/${sample}-ambiguous.1.fq) / 4) >> ${out_dir}/ReadCount/${sample}.txt
# echo $(wc -l < ${out_dir}/${sample}_R1.fq)$'\t'$(wc -l < ${out_dir}/${sample}-graft.1.fq)$'\t'$(wc -l < ${out_dir}/${sample}-host.1.fq)$'\t'$(wc -l < ${out_dir}/${sample}-both.1.fq)$'\t'$(wc -l < ${out_dir}/${sample}-neither.1.fq)$'\t'$(wc -l < ${out_dir}/${sample}-ambiguous.1.fq) >> ${out_dir}/ReadCount/${sample}.txt

# Arrange the outputs from xengsort.
rm ${out_dir}/${sample}_R1.fq
rm ${out_dir}/${sample}_R2.fq

gzip ${out_dir}/${sample}-graft.1.fq
mv ${out_dir}/${sample}-graft.1.fq.gz ${out_dir}/Human/${sample}_R1.fq.gz
ln -s ${out_dir}/Human/${sample}_R1.fq.gz ${human_dir}/${sample}_R1.fq.gz

gzip ${out_dir}/${sample}-graft.2.fq
mv ${out_dir}/${sample}-graft.2.fq.gz ${out_dir}/Human/${sample}_R2.fq.gz
ln -s ${out_dir}/Human/${sample}_R2.fq.gz ${human_dir}/${sample}_R2.fq.gz

gzip ${out_dir}/${sample}-host.1.fq
mv ${out_dir}/${sample}-host.1.fq.gz ${out_dir}/Mouse/${sample}_R1.fq.gz

gzip ${out_dir}/${sample}-host.2.fq
mv ${out_dir}/${sample}-host.2.fq.gz ${out_dir}/Mouse/${sample}_R2.fq.gz

gzip ${out_dir}/${sample}-both.1.fq
mv ${out_dir}/${sample}-both.1.fq.gz ${out_dir}/Both/${sample}_R1.fq.gz

gzip ${out_dir}/${sample}-both.2.fq
mv ${out_dir}/${sample}-both.2.fq.gz ${out_dir}/Both/${sample}_R2.fq.gz

gzip ${out_dir}/${sample}-neither.1.fq
mv ${out_dir}/${sample}-neither.1.fq.gz ${out_dir}/Neither/${sample}_R1.fq.gz

gzip ${out_dir}/${sample}-neither.2.fq
mv ${out_dir}/${sample}-neither.2.fq.gz ${out_dir}/Neither/${sample}_R2.fq.gz

gzip ${out_dir}/${sample}-ambiguous.1.fq
mv ${out_dir}/${sample}-ambiguous.1.fq.gz ${out_dir}/Ambiguous/${sample}_R1.fq.gz

gzip ${out_dir}/${sample}-ambiguous.2.fq
mv ${out_dir}/${sample}-ambiguous.2.fq.gz ${out_dir}/Ambiguous/${sample}_R2.fq.gz

####################################################################################################################
####################################################################################################################
