####################################################################################################################
####################################################################################################################
# Create genotype vcf for each sample with genotype locations decided by dbSNP.
# Author: Haiying Kong
# Last Modified: 10 June 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

# Set parameters.

####################################################################################################################
####################################################################################################################
# Read in common variant locations.
dbSNP = read.table('Reference/FASTQuick/dbSNP/Target1/dbSNP_Target1.vcf', header=FALSE, sep='\t')[ ,1:5]

# Clean or create folder to save genotype vcf files.
geno.dir = 'QC/Lock/VerifyBamID/Genotype/'
if (dir.exists(geno.dir))   {
  unlink(geno.dir, recursive=TRUE)
  }
dir.create(geno.dir, recursive=TRUE)

# Get batch names.
batch.nums = str_pad(1:13, 3, pad='0')
batches = paste0('Primary_', batch.nums)

for (batch in batches)  {
  # Create a folder for this batch to save genotype vcf files.
  dir.create(paste0(geno.dir, '/', batch))
  vcf.dir = paste0('BatchWork/', batch, '/Lock/SNV_InDel/VarDict/vcf/')
  vcf.files = dir(vcf.dir, pattern='.vcf')

  for (vcf.file in vcf.files)  {
    n.rows = system(paste0("grep -v '^#' ", "/home/projects/cu_10184/projects/PTH/", vcf.dir, vcf.file, " | wc -l"), intern=TRUE)
    if (n.rows>0)  {
      vcf = read.table(paste0(vcf.dir, vcf.file), header=FALSE, sep='\t')[ ,c(1,2,4,5,7,8)]
      aster = merge(dbSNP, vcf, by=names(dbSNP)[-3], all.x=TRUE, all.y=FALSE)
      aster$V7[is.na(aster$V7)] = 'PASS'
      aster$V8 = sapply(aster$V8, function(x)    {
                                    if (!is.na(x))  {
                                      tmp = unlist(strsplit(x, ';'))
                                      unlist(tmp[grep('^AF=', tmp)])
                                      }  else  {
                                      'AF=0.0000'
                                      }
                                    }  )
      aster$V9 = '.'
      aster = aster[ ,c(1,2,5,3,4,8,6,7)]
      names(aster) = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
      write.table(aster, paste0(geno.dir, batch, '/', vcf.file), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
      }
    }
  }

####################################################################################################################
####################################################################################################################
