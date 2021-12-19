####################################################################################################################
####################################################################################################################
# Create large table for all variants identified from all samples.
# Author: Haiying Kong
# Last Modified: 19 December 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
# args = commandArgs(trailingOnly=TRUE)
# proj.name = args[1]

library(xlsx)
library(stringr)

# Set parameters.
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

####################################################################################################################
####################################################################################################################
# Create large maf file by concatenating all maf files from all samples.
####################################################################################################################
# Create result folder to save large tables for all samples.
res.dir = '/home/projects/cu_10184/projects/PTH/AllBatches/Result/SNV_InDel/Filtered'
if (dir.exists(res.dir))  unlink(res.dir, recursive=TRUE)
dir.create(res.dir)

for (filter.style in c('Long', 'Medium', 'Short'))  {
  apple = c()

  for (batch in batches)  {
    dir.name = paste0('/home/projects/cu_10184/projects/PTH/BatchWork/', batch, '/Result/SNV_InDel/Filtered/', filter.style)
    maf.files = dir(dir.name, pattern='.maf')

    for (maf.file in maf.files)  {
      aster = read.table(paste(dir.name, maf.file, sep='/'), header=TRUE, quote='', fill=TRUE, sep='\t')
      apple = rbind(apple, aster)
      }
    }

  write.table(apple, paste0(res.dir, filter.style, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  write.xlsx(apple, paste0(res.dir, filter.style, '.xlsx'), sheetName=filter.style,  
             row.names=FALSE, col.names=TRUE, append=FALSE)
  }

####################################################################################################################
####################################################################################################################
