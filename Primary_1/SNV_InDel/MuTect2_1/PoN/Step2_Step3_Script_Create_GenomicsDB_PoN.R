####################################################################################################################
####################################################################################################################
# Create script to create PoN with panel2 normals for Step2 and Step3.
# Author: Haiying Kong
# Last Modified: 22 July 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/PoN')
options(stringsAsFactors=FALSE)
rm(list=ls())

# Step2: Create a GenomicsDB from the normal Mutect2 calls.
aster1 = 'gatk GenomicsDBImport'
aster1 = paste(aster1, '-R /home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta')
aster1 = paste(aster1, '-L /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2.bed')
aster1 = paste(aster1, '--genomicsdb-workspace-path /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2_1/PoN/GenomicsDB')

vcf.dir = '/home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2_1/PoN/Normal'
vcf = dir(vcf.dir, pattern='.vcf.gz$')
vcf = paste0(' -V ', vcf.dir, '/', vcf)
vcf = paste(vcf, collapse='')

aster1 = paste0(aster1, vcf)

apple = rbind('rm -rf /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2_1/PoN/tmp',
              'mkdir /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2_1/PoN/tmp',
              'cd /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2_1/PoN/tmp',
              'rm -rf /home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2_1/PoN/GenomicsDB',
              aster1)
row.names(apple) = NULL

# Step3: Combine the normal calls using CreateSomaticPanelOfNormals.
aster2 = 'gatk CreateSomaticPanelOfNormals'
aster2 = paste(aster2, '-R /home/projects/cu_10184/people/haikon/Reference/GATK/hg38_MaskedU2AF1L5/Homo_sapiens_assembly38_MaskedU2AF1L5.fasta')
aster2 = paste(aster2, '-V gendb:///home/projects/cu_10184/projects/PTH/AllBatches_1/Lock/SNV_InDel/MuTect2_1/PoN/GenomicsDB')
aster2 = paste(aster2, '-O /home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon_1.vcf.gz')
aster2 = paste(aster2, '--min-sample-count 2')

# Concatenate Step2 and Step3.
apple = rbind(apple,
              'rm -f /home/projects/cu_10184/projects/PTH/Reference/MuTect2/pon_1.vcf.gz',
              aster2)
row.names(apple) = NULL

# Save the list of variants after filtering.
write.table(apple, 'Step2_Step3_Create_GenomicsDB_PoN_job.sh', row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# Make the script excutable.
system('chmod +x /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/PoN/Step2_Step3_Create_GenomicsDB_PoN_job.sh')

####################################################################################################################
####################################################################################################################
