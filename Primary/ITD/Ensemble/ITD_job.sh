####################################################################################################################
####################################################################################################################
# Identify ITDs with Pindel and getITD, combine filter results with help from VarDict - job.
# Author: Haiying Kong and Balthasar Schlotmann
# Last Modified: 2 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

#PBS -l nodes=1:ppn=8
#PBS -l mem=20GB
#PBS -l walltime=70:00:00

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate

####################################################################################################################
####################################################################################################################
# Reference databases:
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
target_flt3="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/FLT3.bed"
getITD_reference="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon.txt"
getITD_anno="/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv"

# Software tools:
ScanITD="/home/projects/cu_10184/people/haikon/Software/ScanITD/ScanITD.py"
getITD="python3 /home/projects/cu_10184/people/haikon/Software/getITD/getitd.py"

# Set parameters.
scheme="Scheme_2"
thresh_n_alt=1

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

####################################################################################################################
# VarDict:
####################################################################################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${scheme}/VarDict_FilterClean.R ${batch} ${sample} ${Lock_VarDict_dir} ${Lock_ITD_dir}/VarDict

####################################################################################################################
# Pindel:
####################################################################################################################
# Run Pindel.
conda activate
echo -e ${BAM_dir}/${sample}.bam"\t"250"\t"${sample} > ${Lock_ITD_dir}/Pindel/${sample}_config.txt
pindel -f $hg -t ${target_itd} -c chr13:28033736-28034557  \
       -i ${Lock_ITD_dir}/Pindel/${sample}_config.txt -o ${Lock_ITD_dir}/Pindel/${sample}  \
       -T $n_thread -x 3 -u 0.04
rm ${Lock_ITD_dir}/Pindel/${sample}_config.txt
pindel2vcf -p ${Lock_ITD_dir}/Pindel/${sample}_TD -r $hg -R HG38 -d 20201224 -v ${Lock_ITD_dir}/Pindel/${sample}.vcf
conda deactivate

# Clean the output vcf.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${scheme}/Pindel_Clean.R ${batch} ${sample} ${Lock_ITD_dir}/Pindel

####################################################################################################################
# ScanITD:
####################################################################################################################
# Run ScanITD.
conda activate ScanITD
$ScanITD -r $hg -t ${target_itd} -i ${BAM_dir}/${sample}.bam -o ${Lock_ITD_dir}/ScanITD/${sample} -m 20 -c 1 -d 100 -f 0.0001 -l 5 -n 3
conda deactivate

#  -m MAPQ, --mapq MAPQ  minimal MAPQ in BAM for calling ITD (default: 15)
#  -c AO, --ao AO        minimal observation count for ITD (default: 4)
#  -d DP, --depth DP     minimal depth to call ITD (default: 10)
#  -f VAF, --vaf VAF     minimal variant allele frequency (default: 0.1)
#  -l ITD_LEN, --len ITD_LEN
#                        minimal ITD length to report (default: 10)
#  -n MISMATCH           maximum allowed mismatch bases of pairwise local alignment (default: 3)

# Clean the output tsv.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${scheme}/ScanITD_Clean.R ${batch} ${sample} ${Lock_ITD_dir}/ScanITD

####################################################################################################################
# getITD:
####################################################################################################################
# Filter bam file with target bed file and conver it to fastq.
mkdir ${Lock_ITD_dir}/getITD/${sample}
cd ${Lock_ITD_dir}/getITD/${sample}
bedtools intersect -abam ${BAM_dir}/${sample}.bam -b ${target_flt3} -u > FLT3.bam
# samtools view -b ${BAM_dir}/${sample}.bam chr13:28033736-28034557 > FLT3.bam
samtools sort -n FLT3.bam > tmp.bam
mv tmp.bam FLT3.bam
bedtools bamtofastq -i FLT3.bam -fq FLT3_R1.fq -fq2 FLT3_R2.fq

# Run getITD.
conda activate getITD
$getITD -reference ${getITD_reference} -anno ${getITD_anno} -infer_sense_from_alignment True -plot_coverage True -require_indel_free_primers False -min_read_length 50 -min_bqs 20 -min_read_copies 1 -filter_ins_vaf 0.001 ${sample} FLT3_R1.fq FLT3_R2.fq
conda deactivate

# Clean the output folder.
rm FLT3*
mv ${sample}_getitd/* ./
rm -r ${sample}_getitd

# Clean the output tsv.
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${scheme}/getITD_Clean.R ${batch} ${sample} ${Lock_ITD_dir}/getITD

####################################################################################################################
# IGV:
####################################################################################################################

ln -s ${BAM_dir}/${sample}.bai ${BAM_dir}/${sample}.bam.bai

cat > ${Lock_ITD_dir}/IGV/${sample}_batch.bat << EOF
new 
genome hg38 
snapshotDirectory ${Lock_ITD_dir}/IGV 
load ${BAM_dir}/${sample}.bam
maxPanelHeight -1 
preference SAM.SHOW_SOFT_CLIPPED TRUE
goto chr13:28033736-28034557
snapshot ${sample}.png 
exit
EOF

cat ${Lock_ITD_dir}/IGV/${sample}_batch.bat

/home/projects/cu_10184/people/haikon/Software/IGV_Linux_2.9.0/igv_auto.sh -b ${Lock_ITD_dir}/IGV/${sample}_batch.bat

rm ${BAM_dir}/${sample}.bam.bai
rm ${Lock_ITD_dir}/IGV/${sample}_batch.bat

####################################################################################################################
# Combine the outcome from different tools to come up with final ITD list.
####################################################################################################################
Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${scheme}/CombineFilter_Pindel_getITD_VarDictAnno.R ${batch} ${sample} ${Lock_ITD_dir} ${Result_ITD_dir} ${thresh_n_alt}

####################################################################################################################
####################################################################################################################
