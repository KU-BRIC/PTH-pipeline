####################################################################################################################
####################################################################################################################
# Fill missing DPs and AFs in maf files output from Funcotator.
# Author: Haiying Kong
# Last Modified: 30 March 2021
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

# Read in maf file.
maf = read.table(paste0(maf.dir, sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')

# Find variants with missing DP or AF.
idx = which(is.na(maf$DP) | is.na(maf$AF))

if (length(idx)>0)  {
  # Find DP and AF for these variants from vcf file.
  vcf = read.table(paste0(vcf.dir, sam, '.vcf'), header=FALSE, sep='\t')
  if (!identical(vcf$V1, maf$Chromosome))
    print(paste0(sam, ' ', call.method, ': VCF and MAF do not match'))  else {
    aster = vcf[idx, ]
    ans = sapply(aster$V8, function(x)  {
                             y = unlist(strsplit(x, ';'))
                             dp = as.numeric(sub('DP=', '', y[grep('DP=', y)]))
                             af = as.numeric(sub('AF=', '', y[grep('AF=', y)]))
                             list(c(dp, af))
                             } )
    maf$DP[idx] = as.numeric(sapply(ans, function(x) unlist(x)[1]))
    maf$AF[idx] = as.numeric(sapply(ans, function(x) unlist(x)[2]))
    write.table(maf, paste0(maf.dir, sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
    }
  }

####################################################################################################################
####################################################################################################################
