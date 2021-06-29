####################################################################################################################
####################################################################################################################
# Collet all ITDs from VarDict, Pindel, ScanITD, getITD, and getITD_1.
# Author: Haiying Kong
# Last Modified: 3 June 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)
library(xlsx)

####################################################################################################################
# Set parameter.
thresh.n.alt = 1

####################################################################################################################
####################################################################################################################
# Collect all ITD calls from all samples.
####################################################################################################################
# Get batch list.
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))
  
# Collect ITDs from all samples.
apple = c()
for (batch in batches)  {
  itd.files = dir(paste0('BatchWork/', batch, '/Result/ITD/Table'), pattern='.txt')
  for (itd.file in itd.files)  { 
    itd = read.table(paste0('BatchWork/', batch, '/Result/ITD/Table/', itd.file), header=TRUE, sep='\t')
    if (length(itd)>0)  {
      apple = rbind(apple, itd)
      }
    }
  }



####################################################################################################################
# Save the results.
write.table(apple, 'Ensemble/Final_ITD_List.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.xlsx(apple, 'Ensemble/Final_ITD_List.xlsx', sheetName='getITD',
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
