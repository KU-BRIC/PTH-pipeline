####################################################################################################################
####################################################################################################################
# Create a summary table for read counts.
# Author: Haiying Kong
# Last Modified: 9 March 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Set values.
batch = 'PDX_001'
classify.tool = 'xengsort'

# Get directories.
in.dir = paste0('PDXseqData/', batch, '/', classify.tool, '/ReadCount/')
out.dir = paste0('PanelSeqData/', batch, '/meta/')
if (!dir.exists(out.dir))  dir.create(out.dir)

####################################################################################################################
####################################################################################################################
# Get list of files.
file.names = dir(in.dir, pattern='.txt')

# Concatenate read counts for all samples.
apple = c()
samples = c()
for (file.name in file.names)  {
  one = read.table(paste0(in.dir, file.name), header=FALSE, sep='\t')
  sam = sub('.txt', '', file.name)
  samples = c(samples, sam)
  apple = rbind(apple, one)
  }
names(apple) = c('All', 'Human', 'Mouse', 'Both', 'Neither', 'Ambiguous')
apple$Human_pct = round(apple$Human/apple$All, 4)
apple = cbind(samples, apple)
names(apple)[1] = 'Sample'
apple = apple[order(-apple$Human_pct), ]

write.table(apple, paste0(out.dir, 'Summary_ReadCounts_', classify.tool, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
