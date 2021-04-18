####################################################################################################################
####################################################################################################################
# Get ITDs from all samples.
# Author: Haiying Kong
# Last Modified: 17 April 2021
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
  dir.name = paste0('BatchWork/', batch, '/Lock/ITD/Pindel/')
  file.names = dir(dir.name, pattern='.vcf')

  for (file.name in file.names)  {
    n.rows = system(paste0("grep -v '^#' ", dir.name, file.name, " | wc -l"), intern=TRUE)
    if (n.rows > 0)  {
      aster = read.table(paste0(dir.name, file.name), header=FALSE, sep='\t')[ ,c(1,2,8,10)]
      apple = rbind(apple, cbind(batch, sub('.vcf', '', file.name), aster))
      }
    }
  }
names(apple) = c('Batch', 'Sample', 'Chrom', 'Start', 'Anno', 'Anno1')
write.table(apple, 'AllBatches/ITD/Pindel/Pindel.vcf', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# Clean and sort the result.
aster = sapply(apple[ ,5], function(x)  {
                             tmp = unlist(strsplit(x, ';'))
                             list(tmp)
                             } )
apple$End = sapply(aster, function(x) as.numeric(sub('END=', '', x[1])))
apple$SVLEN = sapply(aster, function(x) as.numeric(sub('SVLEN=', '', x[3])))
apple$HOMLEN = sapply(aster, function(x) as.numeric(sub('HOMLEN=', '', x[2])))
apple$NTLEN = sapply(aster, function(x) as.numeric(sub('NTLEN=', '', x[5])))

aster = sapply(apple[ ,6], function(x)  {
                             tmp = unlist(strsplit(x, ':'))[2]
                             tmp = unlist(strsplit(tmp, ','))
                             list(tmp)
                             } )
apple$N_Ref = sapply(aster, function(x) as.numeric(x[1]))
apple$N_Alt = sapply(aster, function(x) as.numeric(x[2]))

apple = apple[ ,-(5:6)]
write.table(apple, 'AllBatches/ITD/Pindel/Pindel.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
# write.xlsx(apple, 'AllBatches/ITD/Pindel/Pindel.xlsx', sheetName='Pindel',
#            row.names=FALSE, col.names=TRUE, append=FALSE)

# Filter for FLT3 exon 13-15.
apple = apple[apple$Chrom=='chr13', ]
apple = apple[((apple$Start>28033736 & apple$End<28034141) | (apple$Start>28033931 & apple$End<28034364) | (apple$Start>28034150 & apple$End<28034557)), ]

# Label validated ITDs.
apple$Validated = 0
valid = read.table('/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/Validated.txt', header=TRUE, sep='\t')
apple$Patient = sapply(apple$Sample, function(x) unlist(strsplit(x, '-CMP'))[1])
apple$Validated[apple$Patient %in% valid$Patient] = 1
apple = apple[ ,c('Batch', 'Patient', 'Sample', 'Validated', 'Chrom', 'Start', 'End', 'SVLEN', 'N_Alt', 'N_Ref', 'HOMLEN', 'NTLEN')]
apple = apple[order(-apple$N_Alt, apple$Batch, apple$Patient, apple$Sample, apple$Validated, apple$Start), ]

write.table(apple, 'AllBatches/ITD/Pindel/Pindel_FLT3.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.xlsx(apple, 'AllBatches/ITD/Pindel/Pindel_FLT3.xlsx', sheetName='Pindel',
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
