####################################################################################################################
####################################################################################################################
# Clean the vcf file from Pindel.
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
# Check if the sample has unempty vcf file.
vcf.file = paste0(lock.dir, '/', sam, '.vcf')
n.rows = system(paste0("grep -v '^#' ", vcf.file, " | wc -l"), intern=TRUE)

# Clean unempty vcf and save as txt file.
if (n.rows > 0)  {
  apple = read.table(vcf.file, header=FALSE, sep='\t')[ ,c(1,2,8,10)]
  apple$Batch = batch
  apple$Sample = sam
  apple = apple[ ,c(5,6,1:4)]
  names(apple) = c('Batch', 'Sample', 'Chrom', 'Start', 'Anno', 'Anno1')

  # Clean and sort the result.
  aster = sapply(apple[ ,5], function(x)  {
                               tmp = unlist(strsplit(x, ';'))
                               list(tmp)
                               } )
  apple$End = sapply(aster, function(x) as.numeric(sub('END=', '', x[1])))
  apple$SVLEN = sapply(aster, function(x) as.numeric(sub('SVLEN=', '', x[3])))
  apple$HOMLEN = sapply(aster, function(x) as.numeric(sub('HOMLEN=', '', x[2])))
  apple$NTLEN = sapply(aster, function(x) as.numeric(sub('NTLEN=', '', x[5])))

  aster = sapply(apple[ ,6], function(x)  {
                               tmp = unlist(strsplit(x, ':'))[2]
                               tmp = unlist(strsplit(tmp, ','))
                               list(tmp)
                               } )
  apple$N_Alt = sapply(aster, function(x) as.numeric(x[2]))
  apple$N_Ref = sapply(aster, function(x) as.numeric(x[1]))

  apple = apple[ ,-(5:6)]

  # Label validated ITDs.
  apple$Validated = 0
  valid = read.table('/home/projects/cu_10184/projects/PTH/Reference/ITD/FLT3_1/Validated.txt', header=TRUE, sep='\t')
  apple$Patient = sapply(apple$Sample, function(x) unlist(strsplit(x, '-CMP'))[1])
  apple$Validated[apple$Patient %in% valid$Patient] = 1
  apple = apple[ ,c('Batch', 'Patient', 'Sample', 'Validated', 'Chrom', 'Start', 'End', 'SVLEN', 'N_Alt', 'N_Ref', 'HOMLEN', 'NTLEN')]
  apple = apple[order(apple$Batch, apple$Patient, apple$Sample, apple$Validated, apple$Start), ]

  idx = which(apple$N_Alt>1)
  apple = rbind(apple[idx, ], apple[-idx, ])

  write.table(apple, paste0(lock.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
