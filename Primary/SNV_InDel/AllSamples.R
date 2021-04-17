####################################################################################################################
####################################################################################################################
# Create large table for all variants identified from all samples.
# Author: Haiying Kong
# Last Modified: 8 April 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
proj.name = args[1]

library(xlsx)

####################################################################################################################
####################################################################################################################
# Create large maf file by concatenating all maf files from all samples.
####################################################################################################################
# Create result folder to save large tables for all samples.
res.dir = paste0('/home/projects/cu_10184/projects/', proj.name, '/BatchWork/Panel_Trial_YK/Result/SNV_InDel/Filtered/AllSamples/')
if (dir.exists(res.dir))  unlink(res.dir, recursive=TRUE)
dir.create(res.dir)

for (filter.style in c('Long', 'Medium', 'Short'))  {
  apple = c()
  dir.name = paste0('/home/projects/cu_10184/projects/', proj.name, '/BatchWork/Panel_Trial_YK/Result/SNV_InDel/Filtered/', filter.style)
  maf.files = dir(dir.name, pattern='.maf')

  for (maf.file in maf.files)  {
    aster = read.table(paste(dir.name, maf.file, sep='/'), header=TRUE, quote='', sep='\t')
    apple = rbind(apple, aster)
    }

  write.table(aster, paste0(res.dir, filter.style, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  write.xlsx(aster, paste0(res.dir, filter.style, '.xlsx'), sheetName=filter.style,  
             row.names=FALSE, col.names=TRUE, append=FALSE)
  }

####################################################################################################################
####################################################################################################################
