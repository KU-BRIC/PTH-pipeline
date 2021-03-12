####################################################################################################################
####################################################################################################################
# For each sample, concatenate fastq files from different lanes.
# Author: Haiying Kong
# Last Modified: 9 March 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PDX/PanelSeqData')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Set values.
batch = 'PDX_001'
raw.dir = paste0('/home/projects/cu_10184/projects/PTH/PDXseqData/', batch, '/fastq_raw')
new.dir = paste0('/home/projects/cu_10184/projects/PTH/PDXseqData/', batch, '/fastq')

# Clean new.dir.
if (dir.exists(new.dir))  unlink(new.dir, recursive=TRUE)
dir.create(new.dir)

####################################################################################################################
####################################################################################################################
# Get list of raw fastq files.
fastq.files = dir(raw.dir)

# Get list of samples.
samples = sapply(fastq.files, function(x)  unlist(strsplit(x, '_L00'))[1])
samples = unique(samples)

# Concatenate fastq files from different lanes.
for (sam in samples)  {
  sam.fastq.files = fastq.files[grep(sam, fastq.files)]

  r1.files = sam.fastq.files[grep('_R1_001.fastq.gz', sam.fastq.files)]
  r1.files = sort(paste(raw.dir, r1.files, sep='/'))
  system(paste0('cat ', paste(r1.files, collapse=' '), ' > ', new.dir, '/', sam, '_R1.fq.gz'))

  r2.files = sam.fastq.files[grep('_R2_001.fastq.gz', sam.fastq.files)]
  r2.files = sort(paste(raw.dir, r2.files, sep='/'))
  system(paste0('cat ', paste(r2.files, collapse=' '), ' > ', new.dir, '/', sam, '_R2.fq.gz'))
  }

####################################################################################################################
####################################################################################################################
