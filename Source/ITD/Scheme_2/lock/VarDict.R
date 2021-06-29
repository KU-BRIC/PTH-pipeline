####################################################################################################################
####################################################################################################################
# Get list of insertions in chr13:28033736-28034557 identified by VarDict and clean them.
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
vardict.dir = args[3]
lock.dir = args[4]

####################################################################################################################
####################################################################################################################
# Read in maf column names.
maf.cols = read.table('/home/projects/cu_10184/projects/PTH/Reference/MAF_Columns/MAF_Cols_trim0.txt', header=FALSE, sep='\t')[ ,1]
maf.cols = c(maf.cols, 'TYPE', 'SVLEN')

# Read in the list of variants from VarDict call and filter.
maf.file = paste0(vardict.dir, '/maf/', sam, '.maf')
maf = read.table(maf.file, header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols[-c(1:3,12:14,36:114)]]
maf = maf[(maf$Chrom=='chr13' & maf$Start>28033736 & maf$End<28034557), ]
maf = maf[-(which(maf$Variant_Type=='SNP' | maf$Variant_Type=='DEL')), ]

# If any variants left after filtering, clean them and save.
if (nrow(maf)>0)  {
  maf$Batch = batch
  maf$Sample = sam
  maf = maf[ ,c(32,33,2:4,8,15,10,9,30,31)]
  allele2 = maf$Tumor_Seq_Allele2
  maf$Tumor_Seq_Allele2 = 0
  idx = which(maf$TYPE == 'Insertion')
  if (length(idx)>0)  maf$Tumor_Seq_Allele2[idx] = nchar(allele2[idx])
  idx = which(maf$TYPE == 'DUP')
  if (length(idx)>0)  maf$Tumor_Seq_Allele2[idx] = maf$SVLEN[idx]
  maf = maf[ ,-ncol(maf)]
  names(maf) = c('Batch', 'Sample', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', 'Type')
  maf$End = maf$Start + maf$Length - 1
  maf = maf[order(maf$Batch,maf$Sample,maf$Chrom,maf$Start,-maf$N_Alt), ]

  # Save the result.
  write.table(maf, paste0(lock.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
