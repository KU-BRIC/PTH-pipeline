####################################################################################################################
####################################################################################################################
# Get ITDs from all samples.
# Author: Haiying Kong
# Last Modified: 18 April 2021
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
batch.nums = str_pad(1:12, 3, pad='0')
batches = paste0('Primary_', batch.nums)

# Pickup ITDs from all samples.
apple = c()
for (batch in batches)  {
  dir.name = paste0('BatchWork/', batch, '/Lock/ITD/ScanITD/')
  file.names = dir(dir.name, pattern='.vcf')

  for (file.name in file.names)  {
    n.rows = system(paste0("grep -v '^#' ", dir.name, file.name, " | wc -l"), intern=TRUE)
    if (n.rows > 0)  {
      aster = read.table(paste0(dir.name, file.name), header=FALSE, sep='\t')[ ,c(1,2,8)]
      apple = rbind(apple, cbind(batch, sub('.itd.vcf', '', file.name), aster))
      }
    }
  }
names(apple) = c('Batch', 'Sample', 'Chrom', 'Start', 'Anno')

write.table(apple, 'AllBatches/ITD/ScanITD/ScanITD.vcf', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# Clean and sort the result.
aster = sapply(apple[ ,5], function(x)  {
                             tmp = unlist(strsplit(x, ';'))
                             list(tmp)
                             } )
apple$End = sapply(aster, function(x) as.numeric(sub('END=', '', x[9])))
apple$SVLEN = sapply(aster, function(x) as.numeric(sub('SVLEN=', '', x[5])))
apple$N_Alt = sapply(aster, function(x) as.numeric(sub('AO=', '', x[2])))
apple$DP = sapply(aster, function(x) as.numeric(sub('DP=', '', x[3])))
apple$AB = sapply(aster, function(x) as.numeric(sub('AB=', '', x[4])))
apple = apple[ ,-5]

# Aggregate for duplicated rows.
apple = aggregate(cbind(N_Alt,DP,AB) ~ Batch+Sample+Chrom+Start+End+SVLEN, data=apple, FUN=function(x) paste(x,collapse=','))

# Label validated ITDs.
apple$Validated = 0
valid = read.table('/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/Validated.txt', header=TRUE, sep='\t')
apple$Patient = sapply(apple$Sample, function(x) unlist(strsplit(x, '-CMP'))[1])
apple$Validated[apple$Patient %in% valid$Patient] = 1
apple = apple[ ,c('Batch', 'Patient', 'Sample', 'Validated', 'Chrom', 'Start', 'End', 'SVLEN', 'N_Alt', 'DP', 'AB')]
apple = apple[order(apple$Batch, apple$Patient, apple$Sample, apple$Validated, apple$Start), ]

idx = which(apple$N_Alt=='1')
apple = rbind(apple[-idx, ], apple[idx, ])

write.table(apple, 'AllBatches/ITD/ScanITD/ScanITD.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.xlsx(apple, 'AllBatches/ITD/ScanITD/ScanITD.xlsx', sheetName='ScanITD',
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
