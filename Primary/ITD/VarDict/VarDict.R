####################################################################################################################
####################################################################################################################
# Get list of insertions in chr13:28033736-28034557 identified by VarDict.
# Author: Haiying Kong
# Last Modified: 4 May 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)
library(xlsx)

####################################################################################################################
####################################################################################################################
# Set values. 
batch.nums = str_pad(1:12, 3, pad='0')
batches = paste0('Primary_', batch.nums)

# Read in maf column names.
maf.cols = read.table('Reference/MAF_Columns/MAF_Cols_trim0.txt', header=FALSE, sep='\t')[ ,1]
maf.cols = c(maf.cols, 'TYPE')

# Get full list of variants in the region of interest.
apple = c()
for (batch in batches)  {
  # Get directory for the batch.
  maf.dir = paste0('BatchWork/', batch, '/Lock/SNV_InDel/VarDict/maf/')
  
  # Get sample list.
  maf.files = dir(maf.dir, pattern='.maf$')
  samples = gsub('.maf', '', maf.files)
  
  for (sam in samples)  {
    ####################################################
    # Read in the list of variants and filter.
    maf = read.table(paste0(maf.dir, sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols[-c(1:3,12:14)]]
    maf = maf[(maf$Chrom=='chr13' & maf$Start>28033736 & maf$End<28034557), ]
    maf = maf[-(which(maf$Variant_Type=='SNP' | maf$Variant_Type=='DEL')), ]
    if (nrow(maf)>0)  {
      maf = cbind(batch, sam, maf)
      names(maf)[1:2] = c('Batch', 'Sample')
      apple = rbind(apple, maf)
      }
    }
  }

# Save the results.
write.table(apple, 'AllBatches/ITD/VarDict/VarDict.maf', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.xlsx(apple, 'AllBatches/ITD/VarDict/VarDict.xlsx', sheetName='VarDict',
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
