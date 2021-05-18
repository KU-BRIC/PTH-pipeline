####################################################################################################################
####################################################################################################################
# Clean the tsv file from getITD.
# Author: Haiying Kong
# Last Modified: 18 May 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript
  
# Set options and clean the space.
options(stringsAsFactors=FALSE)
rm(list=ls())
    
# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
batch = args[1]
sam = args[2]
lock.dir = args[3]

####################################################################################################################
####################################################################################################################
# Check if the sample has any ITD identified.
sam.dir = paste0(lock.dir, '/', sam)
itd.file = paste0(sam.dir, '/itds_collapsed-is-same_is-similar_is-close_is-same_trailing_hc.tsv')
if (file.exists(itd.file))  {
  apple = read.table(itd.file, header=TRUE, sep='\t')

  # Clean and sort the result.
  apple$Batch = batch
  apple$Chrom = 'chr13'
  cols = c('Batch', 'sample', 'Chrom', 'end_chr13_bp', 'start_chr13_bp', 'domains', 'insertion_site_chr13_bp', 'counts', 'coverage', 'seq')
  apple = apple[ ,cols]
  names(apple) = c('Batch', 'Sample', 'Chrom', 'Start', 'End', 'Domains', 'InsertionSite', 'N_Alt', 'DP', 'Seq')

  # Aggregate for duplicated rows.
  # apple = aggregate(cbind(N_Alt,DP,AB) ~ Batch+Sample+Chrom+Start+End+SVLEN, data=apple, FUN=function(x) paste(x,collapse=','))

  # Label validated ITDs.
  apple$Validated = 0
  valid = read.table('/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/Validated.txt', header=TRUE, sep='\t')
  apple$Patient = sapply(apple$Sample, function(x) unlist(strsplit(x, '-CMP'))[1])
  apple$Validated[apple$Patient %in% valid$Patient] = 1
  apple = apple[ ,c('Batch', 'Patient', 'Sample', 'Validated', 'Chrom', 'Start', 'End', 'Domains', 'InsertionSite', 'N_Alt', 'DP', 'Seq')]
  apple = apple[order(apple$Batch, apple$Patient, apple$Sample, apple$Validated, apple$End), ]

  # Save the results.
  write.table(apple, paste0(lock.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
