####################################################################################################################
####################################################################################################################
# Clean the vcf file from Pindel.
# Author: Haiying Kong
# Last Modified: 25 June 2021
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
vcf.file = paste0(lock.dir, '/', sam, '.vcf')
n.rows = system(paste0("grep -v '^#' ", vcf.file, " | wc -l"), intern=TRUE)

# Clean unempty vcf and save as txt file.
if (n.rows > 0)  {
  apple = read.table(vcf.file, header=FALSE, sep='\t')[ ,c(1,2,8,10)]
  apple$Batch = batch
  apple$Sample = sam
  apple$Hugo_Symbol = 'FLT3'
  apple = apple[ ,c(5:7,1:4)]
  names(apple) = c('Batch', 'Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'Anno', 'Anno1')

  # Clean and sort the result.
  aster = sapply(apple[ ,6], function(x)  {
                               tmp = unlist(strsplit(x, ';'))
                               list(tmp)
                               } )
  apple$End = sapply(aster, function(x) as.numeric(sub('END=', '', x[1])))
  apple$SVLEN = sapply(aster, function(x) as.numeric(sub('SVLEN=', '', x[3])))
  apple$Type = sapply(aster, function(x) sub('SVTYPE=', '', x[4]))
  apple$NTLEN = sapply(aster, function(x) as.numeric(sub('NTLEN=', '', x[5])))
  apple$Length = apple$SVLEN + apple$NTLEN

  aster = sapply(apple[ ,7], function(x)  {
                               tmp = unlist(strsplit(x, ':'))[2]
                               tmp = unlist(strsplit(tmp, ','))
                               list(tmp)
                               } )
  apple$N_Alt = sapply(aster, function(x) as.numeric(x[2]))
  N_Ref = sapply(aster, function(x) as.numeric(x[1]))
  apple$DP = apple$N_Alt + N_Ref
  apple$AF = round(apple$N_Alt/apple$DP, 4)

  apple = apple[ ,-(6:7)]
  apple = apple[ ,c('Batch', 'Sample', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', 'Type', 'SVLEN', 'NTLEN')]
  apple = apple[order(apple$Batch, apple$Sample, apple$Start), ]

  idx = which(apple$N_Alt>1)
  apple = rbind(apple[idx, ], apple[-idx, ])

  write.table(apple, paste0(lock.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
