####################################################################################################################
####################################################################################################################
# Get list of insertions in chr13:28033736-28034557 identified by Soft-Clipping.
# Author: Haiying Kong
# Last Modified: 22 April 2021
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
batch.nums = str_pad(1:11, 3, pad='0')
batches = paste0('Batch', batch.nums)

# Get full list of ITDs.
apple = c()
for (batch in batches)  {
  # Get directory for the batch.
  res.dir = paste0('/home/projects/cu_10145/people/haikon/Project/PTH/BatchWork/', batch, '/Lock/ITD/SoftClipping/')
  
  # Get sample list.
  res.files = dir(res.dir, pattern='.tsv$')
  res.files = res.files[-grep('_spcs.tsv', res.files)]
 
  # Update apple if there is ITD in this batch.
  if (length(res.files) > 0)  { 
    for (res.file in res.files)  {
      aster = read.table(paste0(res.dir, res.file), header=TRUE, quote='', sep='\t')
      apple = rbind(apple, aster)
      }
    }
  }

res.dir = '/home/projects/cu_10184/projects/PTH_test/BatchWork/Panel_Trial_YK/Lock/ITD/SoftClipping/'
res.files = dir(res.dir, pattern='.tsv$')
res.files = res.files[-grep('_spcs.tsv', res.files)]
  
# Update apple if there is ITD in this batch.
if (length(res.files) > 0)  {
  for (res.file in res.files)  {
    aster = read.table(paste0(res.dir, res.file), header=TRUE, quote='', sep='\t')
    apple = rbind(apple, aster)
    }
  }

write.table(apple, 'AllBatches/ITD/SoftClipping/SoftClipping.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.xlsx(apple, 'AllBatches/ITD/SoftClipping/SoftClipping.xlsx', sheetName='SoftClipping',
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
