####################################################################################################################
####################################################################################################################
# Create script to create PoN for Step3.
# Author: Haiying Kong
# Last Modified: 19 July 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2/PoN')
options(stringsAsFactors=FALSE)
rm(list=ls())

# Step3: Combine the normal calls using CreateSomaticPanelOfNormals.
aster2 = 'gatk CreateSomaticPanelOfNormals'
aster2 = paste(aster2, '-R /home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta')
aster2 = paste(aster2, '-V gendb:///home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/GenomicsDB')
aster2 = paste(aster2, '-O /home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon.vcf.gz')

if (file.exists('/home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon.vcf.gz'))
  unlink('/home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon.vcf.gz')

apple = rbind('cd /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2/PoN/tmp', aster2)
row.names(apple) = NULL

# Save the list of variants after filtering.
write.table(apple, 'Step3_Create_PoN.sh', row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# Make the script excutable.
system('chmod +x /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2/PoN/Step3_Create_PoN.sh')

####################################################################################################################
####################################################################################################################
