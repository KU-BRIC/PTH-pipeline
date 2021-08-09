####################################################################################################################
####################################################################################################################
# Get list of variants from all samples.
# Author: Haiying Kong
# Last Modified: 9 August 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH_test')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

####################################################################################################################
####################################################################################################################
# Collect variants from all samples.
####################################################################################################################
####################################################################################################################
# DB:
batch = 'Kirsten_groenbaek_small_panel_data'
scheme.name = 'DB'
long = c()
short = c()

maf.dir = paste0('BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/', scheme.name)

# Long:
maf.files = dir(maf.dir, pattern='_Long.maf')
for (maf.file in maf.files)  {
  maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
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
  maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
  if (nrow(maf)>0)  {
    sam = sub('_Short.maf', '', maf.file)
    pth.id = unlist(strsplit(sam, '-'))[1]
    maf = cbind(batch, pth.id, sam, maf)
    names(maf)[1:3] = c('Batch', 'PTH_ID', 'Sample')
    short = rbind(short, maf)
    }
  }
 
# Save the results.
write.table(long, paste0('BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/AllSamples/', 'FilteredVariants_Long.txt'), 
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

write.table(short, paste0('BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/AllSamples/', 'FilteredVariants_Short.txt'),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
