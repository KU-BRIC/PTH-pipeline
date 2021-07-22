####################################################################################################################
####################################################################################################################
# Get list of Horizon and REDCap variants called by our pipeline.
# Author: Haiying Kong
# Last Modified: 22 July 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)
library(xlsx)


data.name = 'FilteredVariants_VariantClass'

####################################################################################################################
####################################################################################################################
# Read in variant list for all samples.
####################################################################################################################
aster = read.table(paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/', data.name, '.txt'), header=TRUE, quote='', sep='\t')
aster.tag = paste(aster$PTH_ID, aster$Hugo_Symbol, aster$Chromosome, aster$Start_Position, aster$Reference_Allele, aster$Tumor_Seq_Allele2, sep='_')

####################################################################################################################
####################################################################################################################
# Read in Horizon and REDCap variants.
####################################################################################################################
# Horizon:
horizon = read.table('Reference/Horizon/HorizonMyeloid.txt', header=TRUE, sep='\t')[ ,c(1:7,12)]
horizon = cbind(rep('Horizon', nrow(horizon)), horizon)
idx = which(horizon$OnPanel=='yes')
names(horizon)[c(1,8,9)] = c('PTH_ID', 'ProteinChange', 'Note')
horizon$Note[idx] = 'Horizon_OnPanel'
horizon$Note[-idx] = 'Horizon_NotOnPanel'
horizon.tag = paste(horizon$PTH_ID, horizon$Hugo_Symbol, horizon$Chrom, horizon$Start, horizon$Ref, horizon$Alt, sep='_')

# REDCap:
redcap = read.table('Reference/REDCap/REDCap_202102_Anno.txt', header=TRUE, sep='\t')[ ,c(1,2,7,8,10,11,3,12:14,6,5)]
redcap$Note = 'REDCap'
idx = which(redcap$Important==1)
redcap$Note[idx] = paste(redcap$Note[idx], 'Important', sep='_')
idx = which(redcap$SNP==1)
redcap$Note[idx] = paste(redcap$Note[idx], 'SNP', sep='_')
redcap = redcap[ ,-(11:12)]

# Append Horizon and REDCap.
flags = unique(rbind(horizon[ ,1:3], redcap[ ,1:3]))

####################################################################################################################
####################################################################################################################
# Collect all variants that match PTH_ID, Hugo_Symbol and Chrom with Horizon or REDCap.
####################################################################################################################
####################################################################################################################
# Tag with Horizon variants.
####################################################################################################################
idx = match(horizon.tag, aster.tag)
na.idx = which(is.na(idx))

apple = aster[idx[-na.idx], ]
apple$AF_Clinic = horizon$AF[-na.idx]
apple$Note = horizon$Note[-na.idx]


####################################################################################################################
# Tag with REDCap variants.
####################################################################################################################
# Identify REDCap patients that are not in our data.
idx = c(which(!(redcap$PTH_ID %in% unique(aster$PTH_ID))), which(is.na(redcap$cDNA_Change)))
if (length(idx)>0)  {
  write.table(redcap[idx, ], 'AllBatches/SNV_InDel/ClinicTag/REDCap_MissingPatients.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  redcap = redcap[-idx, ]
  }

####################################################################################################################
# Include only REDCap patients in aster.
aster = aster[(aster$PTH_ID %in% redcap$PTH_ID), ]

# Exact match for PTH_ID, Hugo_Symbol, Chrom, Start, Ref and Alt:
redcap.tag = paste(redcap$PTH_ID, redcap$Hugo_Symbol, redcap$Chrom, redcap$Start, redcap$Ref, redcap$Alt, sep='_')
aster.tag = paste(aster$PTH_ID, aster$Hugo_Symbol, aster$Chrom, aster$Start, aster$Ref, aster$Alt, sep='_')
idx = match(redcap.tag, aster.tag)
na.idx = which(is.na(idx))

tmp = aster[idx[-na.idx], ]
tmp$AF_Clinic = redcap$AF[-na.idx]
tmp$Note = redcap$Note[-na.idx]
apple = rbind(apple, tmp)

redcap = redcap[na.idx, ]
aster = aster[-idx[-na.idx], ]

####################################################################################################################
# Exact match for PTH_ID, Hugo_Symbol, Chrom, Start, and cDNA_Change.
redcap.tag = paste(redcap$PTH_ID, redcap$Hugo_Symbol, redcap$Chrom, redcap$Start, redcap$cDNA_Change, sep='_')
aster.tag = paste(aster$PTH_ID, aster$Hugo_Symbol, aster$Chrom, aster$Start, aster$cDNA_Change, sep='_')
idx = match(redcap.tag, aster.tag)
na.idx = which(is.na(idx))

tmp = aster[idx[-na.idx], ]
tmp$AF_Clinic = redcap$AF[-na.idx]
tmp$Note = redcap$Note[-na.idx]
apple = rbind(apple, tmp)

redcap = redcap[na.idx, ]
aster = aster[-idx[-na.idx], ]

####################################################################################################################
# Exact match for PTH_ID, Hugo_Symbol, Chrom, and with loose match for Start.
redcap.start = round(redcap$Start/20)
aster.start = round(aster$Start/20)

redcap.tag = paste(redcap$PTH_ID, redcap$Hugo_Symbol, redcap$Chrom, redcap.start, sep='_')
aster.tag = paste(aster$PTH_ID, aster$Hugo_Symbol, aster$Chrom, aster.start, sep='_')
idx = match(redcap.tag, aster.tag)
na.idx = which(is.na(idx))

redcap.miss = redcap[na.idx, ]

tmp = aster[idx[-na.idx], ]
tmp$AF_Clinic = redcap$AF[-na.idx]
tmp$Note = redcap$Note[-na.idx]
tmp$Start_Clinic = redcap$Start[-na.idx]
tmp$Ref_Clinic = redcap$Ref[-na.idx]
tmp$Alt_Clinic = redcap$Alt[-na.idx]
tmp$Refseq_ID_Clinic = redcap$Refseq_ID[-na.idx]
tmp$cDNA_Change_Clinic = redcap$cDNA_Change[-na.idx]
tmp$Protein_Change_Clinic = redcap$Protein_Change[-na.idx]

tmp = tmp[ ,c(4,1,5,2,3,6,19,7,20,8,21,9,22,10,23,11,24,12:18)]

# Manually check tmp file.
idx = 1:23
tmp = tmp[idx, names(apple)]

apple = rbind(apple, tmp)

apple = apple[ ,c(4,1,5,2,3,6:18)]
apple$AF = round(apple$AF, 3)
apple$AF_VarDict = round(apple$AF_VarDict, 3)
apple$AF_SNVer = round(apple$AF_SNVer, 3)
apple$AF_LoFreq = round(apple$AF_LoFreq, 3)

####################################################################################################################
# Save the results.
write.table(apple, 'AllBatches/SNV_InDel/ClinicTag/CalledReferenceVariants.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.table(redcap.miss, 'AllBatches/SNV_InDel/ClinicTag/MissedREDCapVariants.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

write.xlsx(apple, 'AllBatches/SNV_InDel/ClinicTag/CalledReferenceVariants.xlsx', sheetName='ScanITD',
           row.names=FALSE, col.names=TRUE, append=FALSE)
write.xlsx(redcap.miss, 'AllBatches/SNV_InDel/ClinicTag/MissedREDCapVariants.xlsx', sheetName='ScanITD',
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
