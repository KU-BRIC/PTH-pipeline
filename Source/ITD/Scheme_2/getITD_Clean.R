####################################################################################################################
####################################################################################################################
# Clean the tsv file from getITD.
# Author: Haiying Kong
# Last Modified: 1 June 2021
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
  apple$AF = apple$counts / apple$coverage
  cols = c('Batch', 'sample', 'Chrom', 'end_chr13_bp', 'start_chr13_bp', 'length', 'counts', 'coverage', 'AF', 'trailing')
  apple = apple[ ,cols]
  names(apple) = c('Batch', 'Sample', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', 'Type')
  idx = which(apple$Type == 'True')
  apple$Type = 'non_trailing'
  apple$Type[idx] = 'trailing'

  # Aggregate for duplicated rows.
  # apple = aggregate(cbind(N_Alt,DP,AB) ~ Batch+Sample+Chrom+Start+End+SVLEN, data=apple, FUN=function(x) paste(x,collapse=','))

  # Save the results.
  write.table(apple, paste0(lock.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
