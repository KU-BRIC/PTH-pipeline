####################################################################################################################
####################################################################################################################
# Combine ITD results from Pindel, getITD and update with possible AF and N_Alt from VarDict.
# Author: Haiying Kong
# Last Modified: 22 May 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set working directory, options and clean the space.
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
batch = args[1]
sam = args[2]
lock.dir = args[3]
res.dir = args[4]
thresh.n.alt = args[5]

####################################################################################################################
####################################################################################################################
# Read in ITD calls and organize the data.
####################################################################################################################
# VarDict:
vardict.file = paste0(lock.dir, '/VarDict/', sam, '_clean.maf')
vardict.flag = 0
if (file.exists(vardict.file))  {
  vardict = read.table(vardict.file, header=TRUE, quote='', sep='\t')
  vardict = vardict[vardict$N_Alt>thresh.n.alt, ]
  vardict.flag = 1
  }

####################################################################################################################
# Comnine Pindel and getITD results, then update N_Alt and AF from VarDict.
####################################################################################################################
itd = c()

####################################################################################################################
# Pindel:
pindel.file = paste0(lock.dir, '/Pindel/', sam, '.txt')
if (file.exists(pindel.file))  {
  pindel = read.table(pindel.file, header=TRUE, quote='', sep='\t')
  pindel$Caller = 'Pindel'
  pindel$DP = pindel$N_Alt + pindel$N_Ref
  pindel$AF = round(pindel$N_Alt/pindel$DP, 4)
  pindel$Length = pindel$End - pindel$Start
  pindel = pindel[ ,c('Batch', 'Sample', 'Caller', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', 'Validated')]
  pindel = pindel[pindel$N_Alt>thresh.n.alt, ]
  itd = rbind(itd, pindel)
  }

####################################################################################################################
# getITD:
getitd.file = paste0(lock.dir, '/getITD/', sam, '.txt')
if (file.exists(getitd.file))  {
  getitd = read.table(getitd.file, header=TRUE, quote='', sep='\t')
  getitd$Caller = 'getITD'
  getitd$AF = round(getitd$N_Alt/getitd$DP, 4)
  getitd$Length = getitd$End - getitd$Start
  getitd = getitd[getitd$N_Alt>thresh.n.alt, ]
  getitd = aggregate(cbind(N_Alt,DP,AF) ~ Batch+Sample+Caller+Chrom+Start+End+Length+Validated, data=getitd, max)
  getitd = getitd[ ,c('Batch', 'Sample', 'Caller', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', 'Validated')]
  itd = rbind(itd, getitd)
  }

####################################################################################################################
# Clean combined ITD list from Pindel and getITD.
if (nrow(itd)>0)  {
  itd = itd[order(itd$Batch, itd$Sample, itd$Chrom, itd$End), ]
  itd = aggregate(cbind(Caller,Start,End,Length,N_Alt,DP,AF) ~ Batch+Sample+Chrom+Validated, data=itd, list)

  # Get final list allowing margin by 1~2bp.
  apple = data.frame(Batch='BatchName', Sample='SampleName', Chrom='chr13', Validated=0, Caller='getITD', Start=12345, End=54321, Length=20, N_Alt=12, DP=1234, AF=0.1234)

  for (i in 1:nrow(itd))  {
    caller = unlist(itd$Caller[i])
    start = as.integer(unlist(itd$Start[i]))
    end = as.integer(unlist(itd$End[i]))
    len = as.integer(unlist(itd$Length[i]))
    n_alt = as.integer(unlist(itd$N_Alt[i]))
    dp = as.integer(unlist(itd$DP[i]))
    af = as.numeric(unlist(itd$AF[i]))
    if (length(end)==1)  {
      i.apple = nrow(apple) + 1
      apple[i.apple, ] = NA
      apple[i.apple, 1:4] = itd[i, 1:4]
      apple$Caller[i.apple] = caller
      apple[i.apple, 6:10] = c(start, end, len, n_alt, dp)
      apple$AF[i.apple] = af
      }  else  {
      if (min(diff(end))>2)  {
        nrow.tmp = length(end)
        tmp = data.frame(Batch = rep(itd$Batch[i], nrow.tmp),
                         Sample = rep(itd$Sample[i], nrow.tmp),
                         Chrom = rep(itd$Chrom[i], nrow.tmp),
                         Validated = rep(itd$Validated[i], nrow.tmp),
                         Caller = caller,
                         Start = start,
                         End = end,
                         Length = len,
                         N_Alt = n_alt,
                         DP = dp,
                         AF = af)
        apple = rbind(apple, tmp)
        }  else   {
        i.tmp = 1
        tmp = data.frame(Caller = caller[1],
                         Start = start[1],
                         End = end[1],
                         Length = len[1],
                         N_Alt = n_alt[1],
                         DP = dp[1],
                         AF = af[1])
        for (i.end in 2:length(end))  {
          if ((end[i.end]-end[i.end-1])<=2)  {
            i.tmp = i.tmp + 1
            tmp[i.tmp, ] = NA
            tmp$Caller[i.tmp] = caller[i.end]
            tmp[i.tmp, 2:6] = c(start[i.end], end[i.end], len[i.end], n_alt[i.end], dp[i.end])
            tmp$AF[i.tmp] = af[i.end]
            }  else  {
            tmp.caller = paste(tmp$Caller, collapse=',')
            tmp.start = paste(tmp$Start, collapse=',')
            tmp.end = paste(tmp$End, collapse=',')
            tmp.len = paste(tmp$Length, collapse=',')
            tmp.n_alt = paste(tmp$N_Alt, collapse=',')
            tmp.dp = paste(tmp$DP, collapse=',')
            tmp.af = paste(tmp$AF, collapse=',')
            i.apple = nrow(apple) + 1
            apple[i.apple, ] = NA
            apple[i.apple, 1:4] = itd[i, 1:4]
            apple$Caller[i.apple] = tmp.caller
            apple[i.apple, 6:10] = c(tmp.start, tmp.end, tmp.len, tmp.n_alt, tmp.dp)
            apple$AF[i.apple] = tmp.af

            i.tmp = 1
            tmp = data.frame(Caller = caller[i.end],
                             Start = start[i.end],
                             End = end[i.end],
                             Length = len[i.end],
                             N_Alt = n_alt[i.end],
                             DP = dp[i.end],
                             AF = af[i.end])
            }
          if (i.end == length(end))  {
            if (nrow(tmp)==1)  {
              i.apple = nrow(apple) + 1
              apple[i.apple, ] = NA
              apple[i.apple, 1:4] = itd[i, 1:4]
              apple[i.apple, 5:11] = tmp
              }  else  {
              tmp.caller = paste(tmp$Caller, collapse=',')
              tmp.start = paste(tmp$Start, collapse=',')
              tmp.end = paste(tmp$End, collapse=',')
              tmp.len = paste(tmp$Length, collapse=',')
              tmp.n_alt = paste(tmp$N_Alt, collapse=',')
              tmp.dp = paste(tmp$DP, collapse=',')
              tmp.af = paste(tmp$AF, collapse=',')
              i.apple = nrow(apple) + 1
              apple[i.apple, ] = NA 
              apple[i.apple, 1:4] = itd[i, 1:4]
              apple$Caller[i.apple] = tmp.caller
              apple[i.apple, 6:10] = c(tmp.start, tmp.end, tmp.len, tmp.n_alt, tmp.dp)
              apple$AF[i.apple] = tmp.af
              }
            }
          }
        }
      }
    }
  apple = apple[-1, ]
  for (j in 6:10)  apple[ ,j] = as.character(apple[ ,j])

  ##################################################################################################################
  # Update VAF with VarDict calls if vardict.flag is 1.
  if (vardict.flag==1)  {
    apple$VarDict_N_Alt = 0
    apple$VarDict_DP = 0
    apple$VarDict_AF = 0.0
    apple$VarDict_Seq = ''

    for (i in 1:nrow(apple))  {
      starts = as.numeric(unlist(strsplit(apple$Start[i], ',')))
      var = vardict[(vardict$Batch==apple$Batch[i] & vardict$Sample==apple$Sample[i]), ]
      if (nrow(var)>0)   {
        tmp.n_alt = c()
        tmp.dp = c()
        tmp.af = c()
        tmp.seq = c()
        for (start in starts)  {
          dist = var$Start - start
          if (min(abs(dist)) <= 2)  {
            i.var = which(abs(dist) == min(abs(dist)))
            tmp.n_alt = c(tmp.n_alt, var$N_Alt[i.var])
            tmp.dp = c(tmp.dp, var$DP[i.var])
            tmp.af = c(tmp.af, var$AF[i.var])
            tmp.seq = c(tmp.seq, var$Seq[i.var])
            }
          }
        if (length(tmp.n_alt)>0)   {
          apple$VarDict_N_Alt[i] = paste(unique(tmp.n_alt), collapse=',')
          apple$VarDict_DP[i] = paste(unique(tmp.dp), collapse=',')
          apple$VarDict_AF[i] = paste(unique(tmp.af), collapse=',')
          apple$VarDict_Seq[i] = paste(unique(tmp.seq), collapse=',')
          }
        }
      }
    }

  ##################################################################################################################
  # Save the results.
  apple = apple[order(apple$Batch, apple$Sample, apple$Start), ]
  write.table(apple, paste0(res.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
