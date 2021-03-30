####################################################################################################################
####################################################################################################################
# Copy soft-clipping results from cu_10145 to cu_10184.
# Author: Haiying Kong
# Last Modified: 29 March 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH_test/BatchWork/AllBatches/ITD/SoftClipping')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

####################################################################################################################
# Set values.
batch.nums = str_pad(1:11, 3, pad='0')
batches = paste0('Batch', batch.nums)

####################################################################################################################
# Copy.
for (batch in batches)  {
  dir.name = paste0('/home/projects/cu_10145/people/haikon/Project/PTH/BatchWork/', batch, '/Lock/ITD/SoftClipping')
  file.names = dir(dir.name)
  idx = grep('_spcs.tsv', file.names)
  system(paste0('rm -rf /home/projects/cu_10184/projects/PTH_test/BatchWork/AllBatches/ITD/SoftClipping/', batch))
  system(paste0('mkdir -p /home/projects/cu_10184/projects/PTH_test/BatchWork/AllBatches/ITD/SoftClipping/', batch, '/spcs'))
  for (file.name in file.names)  {
    i = grep('_spcs.tsv', file.name)
    if (length(i)>0)  {
      system(paste0('cp /home/projects/cu_10145/people/haikon/Project/PTH/BatchWork/', batch, '/Lock/ITD/SoftClipping/', file.name, ' /home/projects/cu_10184/projects/PTH_test/BatchWork/AllBatches/ITD/SoftClipping/', batch, '/spcs/'))
      }  else
      system(paste0('cp /home/projects/cu_10145/people/haikon/Project/PTH/BatchWork/', batch, '/Lock/ITD/SoftClipping/', file.name, ' /home/projects/cu_10184/projects/PTH_test/BatchWork/AllBatches/ITD/SoftClipping/', batch, '/'))
    }
  }

####################################################################################################################
####################################################################################################################
