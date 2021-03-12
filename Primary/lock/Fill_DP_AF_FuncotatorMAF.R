####################################################################################################################
####################################################################################################################
# Fill missing DPs and AFs in maf files output from Funcotator.
# Author: Haiying Kong
# Last Modified: 1 March 2021
####################################################################################################################
####################################################################################################################
options(stringsAsFactors=FALSE)
rm(list=ls())

library(prodlim)
library(parallel)

####################################################################################################################
# Set values.
call.methods = c('VarDict', 'SNVer', 'LoFreq')

####################################################################################################################
# Perform parallel computation to fill missing DP and AF in Funcotator output.
for (call.method in call.methods)  {
  # Get the directories.
  vcf.dir = paste0(call.method, '/vcf/')
  maf.dir = paste0(call.method, '/maf/')

  # Get list of samples.
  maf.files = dir(maf.dir, pattern='.maf')
  if (call.method=='SNVer')
    samples = gsub('.snv.maf', '', maf.files[grep('.snv.maf', maf.files)])  else
    samples = gsub('.maf', '', maf.files)

  n.cores = detectCores()
  mclapply(1:length(samples), function(i)  {
    sam = samples[i]

    # For VarDict and LoFreq maf outputs:
    if (call.method != 'SNVer')  {
      # Read in maf file.
      maf = read.table(paste0(maf.dir, sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')

      # Find variants missing DP or AF.
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
      }

    # For SNVer:
    if (call.method == 'SNVer')  {
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
      }
    },
    mc.cores = n.cores)
  }

####################################################################################################################
####################################################################################################################
