####################################################################################################################
####################################################################################################################
# Create script to create GenomicsDB for step2.
# Author: Haiying Kong
# Last Modified: 19 July 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2/PoN')
options(stringsAsFactors=FALSE)
rm(list=ls())

# Clean working directory and change to working directory.
if (dir.exists('/home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/tmp'))  {
  unlink('/home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/tmp', recursive=TRUE)
  dir.create('/home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/tmp')
  }

# Remove the folder where to save output if it exists.
if (dir.exists('/home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/GenomicsDB'))
  unlink('/home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/GenomicsDB', recursive=TRUE)

# Step2: Create a GenomicsDB from the normal Mutect2 calls.
aster1 = 'gatk GenomicsDBImport'
aster1 = paste(aster1, '-R /home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta')
aster1 = paste(aster1, '-L /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2.bed')
aster1 = paste(aster1, '--genomicsdb-workspace-path /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/GenomicsDB')

vcf.dir = '/home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/Normal'
vcf = dir(vcf.dir, pattern='.vcf.gz$')
vcf = paste0(' -V ', vcf.dir, '/', vcf)
vcf = paste(vcf, collapse='')

aster1 = paste0(aster1, vcf)

apple = rbind('cd /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/tmp', aster1)
row.names(apple) = NULL

# Save the list of variants after filtering.
write.table(apple, 'Step2_Create_GenomicsDB.sh', row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# Make the script excutable.
system('chmod +x /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2/PoN/Step2_Create_GenomicsDB.sh')

####################################################################################################################
####################################################################################################################
