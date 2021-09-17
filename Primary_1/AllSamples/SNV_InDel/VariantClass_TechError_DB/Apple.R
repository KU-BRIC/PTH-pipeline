####################################################################################################################
####################################################################################################################
# Create large table for all variants identified from all samples.
# Author: Haiying Kong
# Last Modified: 27 August 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

####################################################################################################################
# Set parameters.
proj.name = 'PTH'
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

####################################################################################################################
# Read in the names of columns that will be kept.
maf.cols = read.table('/home/projects/cu_10184/projects/PTH/Reference/MAF_Columns/MAF_Cols_1.txt', header=FALSE, sep='\t')[ ,1]
maf.cols = c('Batch', 'Sample', maf.cols)

####################################################################################################################
####################################################################################################################
# Create large maf file by concatenating all maf files from all samples.
####################################################################################################################
# Collect all variants from all samples with technical error, variant class and PoN filtered.
long = c()
short = c()

for (batch in batches)  {
  maf.dir = paste0('/home/projects/cu_10184/projects/PTH/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/DB')
  maf.files = dir(maf.dir, pattern='_Long.maf')
  if (length(maf.files)==0)  next
  samples = gsub('_Long.maf', '', maf.files)
  for (sam in samples)  {
    # Long:
    maf = read.table(paste0(maf.dir, '/', sam, '_Long.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')    
    if (nrow(maf)>0)  {
      # Add and trim columns.
      maf = cbind(batch, sam, maf)
      names(maf)[1:2] = c('Batch', 'Sample')
      maf = maf[ ,maf.cols]
      # Update result file.
      long = rbind(long, maf)
      }

    # Short:
    maf = read.table(paste0(maf.dir, '/', sam, '_Short.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')
    if (nrow(maf)>0)  {
      # Add and trim columns.
      maf = cbind(batch, sam, maf)
      names(maf)[1:2] = c('Batch', 'Sample')
      maf = maf[ ,maf.cols]
      # Update result file.
      short = rbind(short, maf)
      }
    }
  }

# Save the results.
write.table(long, '/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/VariantClass_TechError_DB/Long.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.table(short, '/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/Short.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
