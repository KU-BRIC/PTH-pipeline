####################################################################################################################
####################################################################################################################
# Get ITDs from all samples.
# Author: Haiying Kong
# Last Modified: 8 June 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)
library(xlsx)

####################################################################################################################
####################################################################################################################
# Get batch names.
batch.nums = str_pad(1:13, 3, pad='0')
batches = paste0('Primary_', batch.nums)

# Pickup ITDs from all samples.
apple = c()
for (batch in batches[1:10])  {
  dir.name = paste0('BatchWork/', batch, '/Result/ITD/Table/')
  file.names = dir(dir.name, pattern='.txt')

  for (file.name in file.names)  {
    aster = read.table(paste0(dir.name, file.name), header=TRUE, sep='\t')
    if (nrow(aster) > 0)  {
      apple = rbind(apple, aster)
      }
    }
  }

# Save the results.
write.table(apple, 'AllBatches/ITD/AllSamples_AllCallers.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.xlsx(apple, 'AllBatches/ITD/AllSamples_AllCallers.xlsx', sheetName='Pindel',
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
