####################################################################################################################
####################################################################################################################
# Create symbolic links for PDX fastq files to PDX project folder.
# Author: Haiying Kong
# Last Modified: 9 March 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Set values.
batch = 'PDX_001'

original.dir = '/home/projects/cu_10184/data/projects/pth/data/data_raw/fastq/panel_seq/PDX/20210212_FMP007-9_PDX'
new.dir = paste('PDXseqData', batch, 'fastq_raw', sep='/')

####################################################################################################################
####################################################################################################################
# Create new folder to soft link the fastq files.
system(paste0('rm -rf ', new.dir))
system(paste0('mkdir -p ', new.dir))

# Get list of folders for fastq files.
dir.names = dir(original.dir)

# Crete symbolic links for the fastq files in each folder.
for (dir.name in dir.names)  {
  fastq.dir = paste0(original.dir, '/', dir.name)
  fastq.files = dir(fastq.dir, pattern='.fastq.gz')
  for (fastq.file in fastq.files)  {
    system(paste0('ln -s ', fastq.dir, '/', fastq.file, ' ', new.dir))
    }
  }

####################################################################################################################
####################################################################################################################
