####################################################################################################################
####################################################################################################################
# Get list of variants from all samples.
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

####################################################################################################################
####################################################################################################################
# Get batch list.
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

# Read in the names of columns that will be kept.
maf.cols = read.table('/home/projects/cu_10184/projects/PTH/Reference/MAF_Columns/MAF_Cols_full_1.txt', header=FALSE, sep='\t')[ ,1]

####################################################################################################################
####################################################################################################################
# Collect variants from all samples.
####################################################################################################################
####################################################################################################################
# VariantClass:
scheme.name = 'VariantClass'
apple = c()
for (batch in batches)  {
  maf.dir = paste0('BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/', scheme.name)
  maf.files = dir(maf.dir, pattern='.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols]
    if (nrow(maf)>0)  {
      sam = sub('.maf', '', maf.file)
      pth.id = unlist(strsplit(sam, '-'))[1]
      maf = cbind(batch, pth.id, sam, maf)
      names(maf)[1:3] = c('Batch', 'PTH_ID', 'Sample')
      apple = rbind(apple, maf)
      }
    }
  }

# Save the results.
write.table(apple, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '.txt'),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
# write.xlsx(apple, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '.xlsx'),
#            row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
# TechError:
scheme.name = 'TechError'
apple = c()
for (batch in batches)  {
  maf.dir = paste0('BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/', scheme.name)
  maf.files = dir(maf.dir, pattern='.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols]
    if (nrow(maf)>0)  {
      sam = sub('.maf', '', maf.file)
      pth.id = unlist(strsplit(sam, '-'))[1]
      maf = cbind(batch, pth.id, sam, maf)
      names(maf)[1:3] = c('Batch', 'PTH_ID', 'Sample')
      apple = rbind(apple, maf)
      }
    }
  } 
  
# Save the results.
write.table(apple, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '.txt'), 
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
# write.xlsx(apple, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '.xlsx'),
#            row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
# DB:
scheme.name = 'DB'
long = c()
short = c()
for (batch in batches)  {
  maf.dir = paste0('BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/', scheme.name)

  # Long:
  maf.files = dir(maf.dir, pattern='_Long.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols]
    if (nrow(maf)>0)  {
      sam = sub('_Long.maf', '', maf.file)
      pth.id = unlist(strsplit(sam, '-'))[1]
      maf = cbind(batch, pth.id, sam, maf)
      names(maf)[1:3] = c('Batch', 'PTH_ID', 'Sample')
      long = rbind(long, maf)
      }
    }

  # Short:
  maf.files = dir(maf.dir, pattern='_Short.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols]
    if (nrow(maf)>0)  {
      sam = sub('_Short.maf', '', maf.file)
      pth.id = unlist(strsplit(sam, '-'))[1]
      maf = cbind(batch, pth.id, sam, maf)
      names(maf)[1:3] = c('Batch', 'PTH_ID', 'Sample')
      short = rbind(short, maf)
      }
    }
  }

 
# Save the results.
write.table(long, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '_Long.txt'), 
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
# write.xlsx(long, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '_Long.xlsx'),
#            row.names=FALSE, col.names=TRUE, append=FALSE)

write.table(short, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '_Short.txt'),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
# write.xlsx(short, paste0('AllBatches_1/Result/SNV_InDel/MuTect2_1/FilteredVariants_', scheme.name, '_Short.xlsx'),
#            row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
