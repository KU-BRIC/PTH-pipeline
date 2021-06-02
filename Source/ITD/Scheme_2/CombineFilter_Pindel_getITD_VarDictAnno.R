####################################################################################################################
####################################################################################################################
# Combine ITD results from VarDict, Pindel, ScanITD, getITD.
# Author: Haiying Kong
# Last Modified: 2 June 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set working directory, options and clean the space.
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
batch = args[1]
sam = args[2]
lock.dir = args[3]
res.dir = args[4]
thresh.n.alt = args[5]

####################################################################################################################
####################################################################################################################
# Read in ITD calls from tools and combine them.
####################################################################################################################

apple = c()

####################################################################################################################
# VarDict:
aster.file = paste0(lock.dir, '/VarDict/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  aster = aster[aster$N_Alt>thresh.n.alt, ]
  aster = cbind('VarDict', aster)
  names(aster)[1] = 'Caller'
  apple = rbind(apple, aster)
  }

####################################################################################################################
# Pindel:
aster.file = paste0(lock.dir, '/Pindel/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  aster = aster[aster$N_Alt>thresh.n.alt, ]
  aster = cbind('Pindel', aster)
  names(aster)[1] = 'Caller'
  apple = rbind(apple, aster)
  }

####################################################################################################################
# ScanITD:
aster.file = paste0(lock.dir, '/ScanITD/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  aster = aster[aster$N_Alt>thresh.n.alt, ]
  aster = cbind('ScanITD', aster)
  names(aster)[1] = 'Caller'
  apple = rbind(apple, aster)
  }

####################################################################################################################
# getITD:
aster.file = paste0(lock.dir, '/getITD/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  aster = aster[aster$N_Alt>thresh.n.alt, ]
  aster = cbind('getITD', aster)
  names(aster)[1] = 'Caller'
  apple = rbind(apple, aster)
  }

####################################################################################################################
# Save the results for the sample if the ITD list is not empty.
if (nrow(apple)>0)  {
  apple = apple[order(apple$Batch, apple$Sample, apple$Start), ]
  write.table(apple, paste0(res.dir, '/Table/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
