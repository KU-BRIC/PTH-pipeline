####################################################################################################################
####################################################################################################################
# Create tables for sample pairs for fingerprinting.
# Author: Haiying Kong
# Last Modified: 10 June 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

dir.name = 'PTH'

####################################################################################################################
####################################################################################################################
# Get batch list.
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

####################################################################################################################
# Collect samples from all batches.
all.sam = c()
for (batch in batches)  {
  fastq.files = dir(paste0('/home/projects/cu_10184/projects/', dir.name, '/PanelSeqData/', batch, '/fastq'), pattern='.fq.gz')
  samples = fastq.files[grep('_R1.fq.gz', fastq.files)]
  samples = gsub('_R1.fq.gz', '', samples)
  patients = sapply(samples, function(x) unlist(strsplit(x, '-'))[1])
  all.sam = rbind(all.sam, cbind(patients, batch, samples))
  }

all.sam = as.data.frame(all.sam)
row.names(all.sam) = NULL
names(all.sam) = c('Patient', 'Batch', 'Sample')

all.sam = all.sam[order(all.sam$Patient, all.sam$Sample), ]

####################################################################################################################
# Create table for matched samples.
duplicated.patients = unique(all.sam$Patient[duplicated(all.sam$Patient)])
acorn = all.sam[(all.sam$Patient %in% duplicated.patients), ]

aster = acorn[!duplicated(acorn$Patient), ]
names(aster) = c('Patient', 'Batch_1', 'Sample_1')
acorn = acorn[duplicated(acorn$Patient), ]
names(acorn) = c('Patient', 'Batch_2', 'Sample_2')

apple = c()
for (i in 1:nrow(aster))  {
  idx = which(acorn$Patient == aster$Patient[i])
  n = length(idx)
  apple = rbind(apple, cbind(rep(aster$Patient[i],n), rep(aster$Batch_1[i],n), rep(aster$Sample_1[i],n), acorn[idx, -1]))
  }

apple = data.frame(apple)
names(apple) = c('Patient', 'Batch_1', 'Sample_1', 'Batch_2', 'Sample_2')

apple = apple[ ,-1]
write.table(apple, '/home/projects/cu_10184/projects/PTH/QC/Meta/MatchedSampleList.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
# Create table for unmatched samples.
acorn = all.sam[-((all.sam$Sample %in% apple$Sample_1) | (all.sam$Sample %in% apple$Sample_2)), ]
idx = sample(1:nrow(acorn), 100)
apple = data.frame(Batch_1 = acorn$Batch[idx[1:50]],
                   Sample_1 = acorn$Sample[idx[1:50]],
                   Batch_2 = acorn$Batch[idx[51:100]],
                   Sample_2 = acorn$Sample[idx[51:100]])
write.table(apple, '/home/projects/cu_10184/projects/PTH/QC/Meta/UnmatchedSampleList.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
