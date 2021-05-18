####################################################################################################################
####################################################################################################################
# Find common variants from our variant calls.
# Author: Haiying Kong
# Last Modified: 7 May 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

# Set parameters.
maf.cols = c('Batch', 'Sample', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Strand',
             'Reference_Allele', 'Tumor_Seq_Allele2', 'Start_vcf', 'Ref_vcf', 'Alt_vcf', 'AF', 'DP')
maf.cols = c('Batch', 'Sample', 'Hugo_Symbol', 'Chromosome', 'Start_vcf', 'Strand', 'Ref_vcf', 'Alt_vcf', 'AF', 'DP')
thresh.dp = 100
thresh.freq.min = 0.55
thresh.freq.max = 0.75

####################################################################################################################
####################################################################################################################
# Get batch names.
batch.nums = str_pad(1:12, 3, pad='0')
batches = paste0('Primary_', batch.nums)

# Pickup ITDs from all samples.
apple = c()
n.sam = 0

for (batch in batches)  {
  maf.dir = paste0('BatchWork/', batch, '/Result/SNV_InDel/AllVariants/Callers_Wide/')
  maf.files = dir(maf.dir)
  n.sam = n.sam + length(maf.files)

  for (maf.file in maf.files)  {
    aster = read.table(paste0(maf.dir, maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
    aster = aster[(aster$Variant_Type=='SNP' & aster$DP>thresh.dp), maf.cols]
    aster = aster[!is.na(aster$Batch), -match(c('AF','DP'), names(aster))]
    aster$Freq = 1
    aster = aggregate(Freq~Hugo_Symbol+Chromosome+Start_vcf+Strand+Ref_vcf+Alt_vcf, data=aster, sum)
    apple = rbind(apple, aster)
    }
  }

apple = aggregate(Freq~., data=apple, sum)

# Filter for common variants.
acorn = apple[(apple$Freq>(n.sam*thresh.freq.min) & apple$Freq<(n.sam*thresh.freq.max)), ]

# Save the results.
write.table(acorn, '/home/projects/cu_10184/projects/PTH/Reference/FASTQuick/CommonVariants.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
