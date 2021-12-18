####################################################################################################################
####################################################################################################################
# Get QC summary table for all sams.
# Author: Haiying Kong
# Last Modified: 18 December 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

####################################################################################################################
####################################################################################################################
# Get batch names.
batch.nums = str_pad(1:13, 3, pad='0')
batches = paste0('Primary_', batch.nums)

# Collect summary tables from all sams.
apple = c()
for (batch in batches)  {
  batch.dir = paste0('/home/projects/cu_10184/projects/PTH/QC/Result/FASTQuick/ByBatch/', batch)
  sams = dir(batch.dir)

  for (sam in sams)  {
    aster = read.table(paste0(batch.dir, '/', sam, '/OurSummary.txt'), header=TRUE, sep='\t')
    aster = cbind(batch, sam, aster)
    names(aster)[1:2] = c('Batch', 'Sample')
    apple = rbind(apple, aster)
    }
  }

# Save the results.
write.table(apple, '/home/projects/cu_10184/projects/PTH/QC/Result/FASTQuick/AllBatches/QC_Summary.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
