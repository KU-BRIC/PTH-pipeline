####################################################################################################################
####################################################################################################################
# Fill missing DPs and AFs in maf files output from Funcotator.
# Author: Haiying Kong
# Last Modified: 29 March 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
lock.dir = args[1]
sam = args[2]
call.method = args[3]

setwd(lock.dir)

####################################################################################################################
# Fill missing DP and AF in Funcotator output.
####################################################################################################################
# Get the directories.
vcf.dir = paste0(lock.dir, '/vcf/')
maf.dir = paste0(lock.dir, '/maf/')
      
# snv:
# Read in maf file.
maf = read.table(paste0(maf.dir, sam, '.snv.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')
idx = which(is.na(maf$DP) | is.na(maf$AC))
if (length(idx)>0)  {
  vcf = read.table(paste0(vcf.dir, sam, '.filter.vcf'), header=FALSE, sep='\t')
  if (!identical(vcf$V1, maf$Chromosome))
    print(paste0(sam, ' ', call.method, ' SNV: VCF and MAF do not match'))  else {
    aster = vcf[idx, ]
    ans = sapply(aster$V8, function(x)  {
                             y = unlist(strsplit(x, ';'))
                             dp = as.numeric(sub('DP=', '', y[grep('DP', y)]))
                             ac = as.numeric(sub('AC=', '', y[grep('AC', y)]))
                             list(c(dp, ac))
                             } )
    maf$DP[idx] = as.numeric(sapply(ans, function(x) unlist(x)[1]))
    maf$AC[idx] = as.numeric(sapply(ans, function(x) unlist(x)[2]))
    write.table(maf, paste0(maf.dir, sam, '.snv.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
    }
  }

# indel:
# Read in maf file.
maf = read.table(paste0(maf.dir, sam, '.indel.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')
idx = which(is.na(maf$DP) | is.na(maf$AC))
if (length(idx)>0)  {
  vcf = read.table(paste0(vcf.dir, sam, '.indel.filter.vcf'), header=FALSE, sep='\t')
  if (!identical(vcf$V1, maf$Chromosome))
    print(paste0(sam, ' ', call.method, ' InDel: VCF and MAF do not match'))  else {
    aster = vcf[idx, ]
    ans = sapply(aster$V8, function(x)  {
                             y = unlist(strsplit(x, ';'))
                             dp = as.numeric(sub('DP=', '', y[grep('DP', y)]))
                             ac = as.numeric(sub('AC=', '', y[grep('AC', y)]))
                             list(c(dp, ac))
                             } )
    maf$DP[idx] = as.numeric(sapply(ans, function(x) unlist(x)[1]))
    maf$AC[idx] = as.numeric(sapply(ans, function(x) unlist(x)[2]))
    write.table(maf, paste0(maf.dir, sam, '.indel.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
    }
  }

####################################################################################################################
####################################################################################################################
