####################################################################################################################
####################################################################################################################
# For each sample, concatenate fastq files from different lanes.
# Author: Haiying Kong
# Last Modified: 21 May 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH/PanelSeqData/Primary_013/fastq')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
####################################################################################################################
# Raw sequence data directory:
raw.dir = "/home/projects/cu_10184/data/projects/pth/data/data_raw/fastq/panel_seq/primary/210423_NB501508_0814_AHMM5WAFX2"

# Create temp folder and softlink raw fastq files.
if (length(dir())>0)  unlink('*', recursive=TRUE)
dir.create('temp')
fastq.dirs = dir(raw.dir)
for (fastq.dir in fastq.dirs)   {
  system(paste0('ln -s ', raw.dir, '/', fastq.dir, '/*.fastq.gz temp'))
  }

# Concatanate the fastq files for each sample.
fastq.files = dir('temp')
fastq.files = fastq.files[grep('PTH', fastq.files)]
samples = sapply(fastq.files, function(x)  unlist(strsplit(x, '_L00'))[1])
samples = unique(samples)

# Concatenate fastq files from different lanes.
for (sam in samples)  {
  new.sam.name = sub('-X01B-', '_', sam)
  new.sam.name = sub('-P001-D01-', '_', sam)

  sam.fastq.files = fastq.files[grep(sam, fastq.files)]

  r1.files = sam.fastq.files[grep('_R1_001.fastq.gz', sam.fastq.files)]
  r1.files = paste0('temp/', sort(r1.files))
  system(paste0('cat ', paste(r1.files, collapse=' '), ' > ', new.sam.name, '_R1.fq.gz'))

  r2.files = sam.fastq.files[grep('_R2_001.fastq.gz', sam.fastq.files)]
  r2.files = paste0('temp/', sort(r2.files))
  system(paste0('cat ', paste(r2.files, collapse=' '), ' > ', new.sam.name, '_R2.fq.gz'))
  }

unlink('temp', recursive=TRUE)

####################################################################################################################
####################################################################################################################
