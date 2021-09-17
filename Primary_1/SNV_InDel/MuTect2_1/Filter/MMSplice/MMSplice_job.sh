####################################################################################################################
####################################################################################################################
# Get MMSplice scores for SNV-InDels in SpliceSite - job.
# Author: Haiying Kong
# Last Modified: 8 September 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

source /home/projects/cu_10184/projects/PTH/Software/envsetup
conda deactivate
conda activate kipoi-gpu-MMSplice__deltaLogitPSI

####################################################################################################################
# Reference files.
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
anno_gtf="/home/projects/cu_10184/people/haikon/Reference/Ensembl/Homo_sapiens.GRCh38.85.chr.uniq_exon.gtf"

####################################################################################################################
####################################################################################################################
# Script for one job:
####################################################################################################################
# Change to working directory.
cd ${temp_dir}

# Result directory:
res_dir=${batch_dir}/Result/SNV_InDel/MuTect2_1/MMSplice

####################################################################################################################
# Prepare vcf file.
vcf_0=${batch_dir}/Lock/SNV_InDel/MuTect2_1/vcf/${sam}.vcf
# Make one line for one variant.
bcftools norm -m-both ${batch_dir}/Lock/SNV_InDel/MuTect2_1/vcf/${sam}.vcf -o ${res_dir}/${sam}_temp1.vcf
# Left-normalization format.
bcftools norm -f ${hg}  ${res_dir}/${sam}_temp1.vcf -o ${res_dir}/${sam}_temp2.vcf

# Run MMSplice.
hg="/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta"
anno_gtf="/home/projects/cu_10184/people/haikon/Reference/Ensembl/Homo_sapiens.GRCh38.85.chr.uniq_exon.gtf.gz"
res_dir=/home/projects/cu_10184/projects/PTH/BatchWork_1/Primary_001/Result/SNV_InDel/MuTect2_1/MMSplice
sam=PTH0008-CMP001

args='{"gtf": "'${anno_gtf}'", "fasta_file": "'$hg'", "vcf_file": "'${res_dir}'/'${sam}'_temp2.vcf", "exon_cut_l": 0, "exon_cut_r": 0, "acceptor_intron_cut": 6, "donor_intron_cut": 6, "acceptor_intron_len": 50, "acceptor_exon_len": 3, "donor_exon_len": 5, "donor_intron_len": 13}'
kipoi predict MMSplice/deltaLogitPSI --dataloader_args=$args -o ${res_dir}/${sam}.tsv


kipoi predict MMSplice/deltaLogitPSI --dataloader_args='{"gtf": "/home/projects/cu_10184/people/haikon/Reference/Ensembl/Homo_sapiens.GRCh38.85_chr.gtf", "fasta_file": "/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta", "vcf_file": "/home/projects/cu_10184/projects/PTH/BatchWork_1/Primary_001/Result/SNV_InDel/MuTect2_1/MMSplice/PTH0008-CMP001_temp2.vcf", "exon_cut_l": 0, "exon_cut_r": 0, "acceptor_intron_cut": 6, "donor_intron_cut": 6, "acceptor_intron_len": 50, "acceptor_exon_len": 3, "donor_exon_len": 5, "donor_intron_len": 13}' -o test.tsv

kipoi predict MMSplice/deltaLogitPSI --dataloader_args='{"gtf": "/home/projects/cu_10184/people/haikon/Reference/Ensembl/Homo_sapiens.GRCh38.85.chr.gtf", "fasta_file": "/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta", "vcf_file": "/home/projects/cu_10184/projects/PTH/BatchWork_1/Primary_001/Result/SNV_InDel/MuTect2_1/MMSplice/PTH0008-CMP001_temp2.vcf", "exon_cut_l": 0, "exon_cut_r": 0, "acceptor_intron_cut": 6, "donor_intron_cut": 6, "acceptor_intron_len": 50, "acceptor_exon_len": 3, "donor_exon_len": 5, "donor_intron_len": 13}' -o test.tsv


kipoi predict MMSplice/deltaLogitPSI \
  --dataloader_args='{"gtf": "gtf_chr", "fasta_file": "fasta_file_chr", "vcf_file": "test.vcf", "exon_cut_l": 0, "exon_cut_r": 0, "acceptor_intron_cut": 6, "donor_intron_cut": 6, "acceptor_intron_len": 50, "acceptor_exon_len": 3, "donor_exon_len": 5, "donor_intron_len": 13}' \
  -o 'test.tsv'

kipoi predict MMSplice/deltaLogitPSI --dataloader_args='{"gtf": "/home/projects/cu_10184/people/haikon/Reference/Ensembl/Homo_sapiens.GRCh37.75.chr.uniq_exon.gtf", "fasta_file": "/home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta", "vcf_file": "test.vcf", "exon_cut_l": 0, "exon_cut_r": 0, "acceptor_intron_cut": 6, "donor_intron_cut": 6, "acceptor_intron_len": 50, "acceptor_exon_len": 3, "donor_exon_len": 5, "donor_intron_len": 13}' -o test.tsv


kipoi predict MMSplice/deltaLogitPSI --dataloader_args='{"gtf": "/home/projects/cu_10184/people/haikon/Reference/Ensembl/Homo_sapiens.GRCh38.85.gtf", "fasta_file": "fasta_file", "vcf_file": "vcf_file", "exon_cut_l": 0, "exon_cut_r": 0, "acceptor_intron_cut": 6, "donor_intron_cut": 6, "acceptor_intron_len": 50, "acceptor_exon_len": 3, "donor_exon_len": 5, "donor_intron_len": 13}' -o test.tsv





####################################################################################################################
####################################################################################################################
