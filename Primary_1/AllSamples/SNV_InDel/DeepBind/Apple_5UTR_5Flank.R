####################################################################################################################
####################################################################################################################
# Create large table for all variants in 5'UTR from all samples.
# Author: Haiying Kong
# Last Modified: 4 September 2021
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
apple = c()

for (batch in batches)  {
  maf.dir = paste0('/home/projects/cu_10184/projects/PTH/BatchWork_1/', batch, '/Lock/SNV_InDel/MuTect2_1/maf')
  maf.files = dir(maf.dir, pattern='.maf')
  if (length(maf.files)==0)  next
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
    maf = maf[(maf$Variant_Classification=="5'Flank" | maf$Variant_Classificatio=="5'UTR"), ]
    if (nrow(maf)==0)  next

    # Add and trim columns.
    maf = cbind(batch, sub('.maf', '', maf.file), maf)
    names(maf)[1:2] = c('Batch', 'Sample')
    maf = maf[ ,maf.cols]

    # Update result file.
    apple = rbind(apple, maf)
    }
  }

# Save the results.
write.table(apple, '/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/DeepBind/UTR5_Flank5.maf',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
