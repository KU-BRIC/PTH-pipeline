####################################################################################################################
####################################################################################################################
# Clean the vcf file from ScanITD.
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
# Check if the sample has unempty vcf file.
vcf.file = paste0(lock.dir, '/', sam, '.itd.vcf')
n.rows = system(paste0("grep -v '^#' ", vcf.file, " | wc -l"), intern=TRUE)

# Clean unempty vcf and save as txt file.
if (n.rows > 0)  {
  apple = read.table(vcf.file, header=FALSE, sep='\t')[ ,c(2,8)]
  apple$Batch = batch
  apple$Sample = sam
  apple$Hugo_Symbol = 'FLT3'
  apple$Chrom = 'chr13'
  apple = apple[ ,c(3:6,1:2)]
  names(apple) = c('Batch', 'Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'Anno')

  # Clean and sort the result.
  aster = sapply(apple[ ,6], function(x)  {
                               tmp = unlist(strsplit(x, ';'))
                               list(tmp)
                               } )
  apple$End = sapply(aster, function(x) as.numeric(sub('END=', '', x[9])))
  apple$Length = sapply(aster, function(x) as.numeric(sub('SVLEN=', '', x[5])))
  apple$N_Alt = sapply(aster, function(x) as.numeric(sub('AO=', '', x[2])))
  apple$DP = sapply(aster, function(x) as.numeric(sub('DP=', '', x[3])))
  apple$AF = sapply(aster, function(x) as.numeric(sub('AB=', '', x[4])))
  apple$Type = sapply(aster, function(x) sub('SVTYPE=', '', x[6]))

  apple = apple[ ,-6]
  apple = apple[ ,c('Batch', 'Sample', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', 'Type')]
  apple = apple[order(apple$Batch, apple$Sample, apple$Start), ]

  idx = which(apple$N_Alt>1)
  apple = rbind(apple[idx, ], apple[-idx, ])

  write.table(apple, paste0(lock.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
