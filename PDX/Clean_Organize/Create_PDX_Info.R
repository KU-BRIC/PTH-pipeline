####################################################################################################################
####################################################################################################################
# Create sample information table for PDX samples.
# Author: Haiying Kong
# Last Modified: 23 October 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# PDX sample information to update:
apple = c()

# Batch name:
batch = 'PDX_001'

# Samples:
sams = dir(paste0('/home/projects/cu_10184/projects/PTH/PanelSeqData/', batch, '/fastq/xengsort'), pattern='_R1.fq.gz')
sams = gsub('_R1.fq.gz', '', sams)

# Mice:
mice = sapply(sams, function(x)  unlist(strsplit(x, '-'))[1])
names(mice) = NULL

# Create sample information for this batch.
aster = data.frame(Batch_ID = batch,
                   Mouse_ID = mice,
                   Sample_ID = sams)

# Update the PDX sample information.
apple = rbind(apple, aster)

# Save the results.
write.table(apple, '/home/projects/cu_10184/projects/PTH/Meta/PDX_Info.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
