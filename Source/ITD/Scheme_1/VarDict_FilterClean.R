####################################################################################################################
####################################################################################################################
# Get list of insertions in chr13:28033736-28034557 identified by VarDict and clean them.
# Author: Haiying Kong
# Last Modified: 22 May 2021
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
maf.cols = c(maf.cols, 'TYPE')

# Read in the list of variants from VarDict call and filter.
maf.file = paste0(vardict.dir, '/maf/', sam, '.maf')
maf = read.table(maf.file, header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols[-c(1:3,12:14)]]
maf = maf[(maf$Chrom=='chr13' & maf$Start>28033736 & maf$End<28034557), ]
maf = maf[-(which(maf$Variant_Type=='SNP' | maf$Variant_Type=='DEL')), ]

# If any variants left after filtering, clean them and save.
if (nrow(maf)>0)  {
  maf$Batch = batch
  maf$Sample = sam
  maf = maf[ ,c(110,111,1:109)]
  write.table(maf, paste0(lock.dir, '/', sam, '_raw.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  maf = maf[ ,c(1,2,4,5,17,12,11,10)]
  names(maf) = c('Batch', 'Sample', 'Chrom', 'Start', 'N_Alt', 'DP', 'AF', 'Seq')
  maf = maf[ ,c('Batch', 'Sample', 'Chrom', 'Start', 'N_Alt', 'DP', 'AF', 'Seq')]
  maf = maf[order(maf$Batch,maf$Sample,maf$Chrom,maf$Start), ]
  maf = aggregate(cbind(Start,N_Alt,DP,AF,Seq) ~ Batch+Sample+Chrom, data=maf, list)
  
  # Aggregate redundancies.
  aster = data.frame(Batch='BatchName', Sample='SampleName', Chrom='chr13', Start=12345, N_Alt=12, DP=1234, AF=0.1234, Seq='GCTA')
  for (i in 1:nrow(maf))  {
    start = as.numeric(unlist(maf$Start[i]))
    n_alt = as.numeric(unlist(maf$N_Alt[i]))
    dp = as.numeric(unlist(maf$DP[i]))
    af = as.numeric(unlist(maf$AF[i]))
    seq = unlist(maf$Seq[i])
    if (length(start)==1)  {
      i.aster = nrow(aster) + 1
      aster[i.aster, ] = NA 
      aster[i.aster, 1:3] = maf[i, 1:3]
      aster[i.aster, 4:6] = c(start, n_alt, dp)
      aster$AF[i.aster] = af
      aster$Seq[i.aster] = seq
      }  else  {
      if (min(diff(start))>2)  {
        nrow.tmp = length(start)
        tmp = data.frame(Batch = rep(maf$Batch[i], nrow.tmp),
                         Sample = rep(maf$Sample[i], nrow.tmp),
                         Chrom = rep(maf$Chrom[i], nrow.tmp),
                         Start = start,
                         N_Alt = n_alt,
                         DP = dp,
                         AF = af,
                         Seq = seq)
        aster = rbind(aster, tmp)
        }  else   {
        i.tmp = 1 
        tmp = data.frame(Start = start[1],
                         N_Alt = n_alt[1],
                         DP = dp[1], 
                         AF = af[1],
                         Seq = seq[1]) 
        for (i.start in 2:length(start))  {
          if ((start[i.start]-start[i.start-1])<=2)  {
            i.tmp = i.tmp + 1
            tmp[i.tmp, ] = NA
            tmp[i.tmp, 1:3] = c(start[i.start], n_alt[i.start], dp[i.start])
            tmp[i.tmp, 4] = af[i.start]
            tmp[i.tmp, 5] = seq[i.start]
            }  else  { 
            tmp.start = median(tmp$Start)
            tmp.n_alt = max(tmp$N_Alt)
            tmp.dp = max(tmp$DP)
            tmp.af = max(tmp$AF)
            tmp.seq = max(tmp$Seq)
            i.aster = nrow(aster) + 1
            aster[i.aster, ] = NA 
            aster[i.aster, 1:3] = maf[i, 1:3]
            aster[i.aster, 4:6] = c(tmp.start, tmp.n_alt, tmp.dp)
            aster$AF[i.aster] = tmp.af
            aster$Seq[i.aster] = tmp.seq
 
            i.tmp = 1
            tmp = data.frame(Start = start[i.start],
                             N_Alt = n_alt[i.start],
                             DP = dp[i.start], 
                             AF = af[i.start],
                             Seq = seq[i.start]) 
            }
          if (i.start == length(start))  {
            if (nrow(tmp)==1)  {
              i.aster = nrow(aster) + 1
              aster[i.aster, ] = NA 
              aster[i.aster, 1:3] = maf[i, 1:3]
              aster[i.aster, 4:8] = tmp
              }  else  {
              tmp.start = median(tmp$Start)
              tmp.n_alt = max(tmp$N_Alt)
              tmp.dp = max(tmp$DP)
              tmp.af = max(tmp$AF)
              tmp.seq = max(tmp$Seq)
              i.aster = nrow(aster) + 1
              aster[i.aster, 1:3] = maf[i, 1:3]
              aster[i.aster, 4:6] = c(tmp.start, tmp.n_alt, tmp.dp)
              aster$AF[i.aster] = tmp.af
              aster$Seq[i.aster] = tmp.seq
              }
            } 
          } 
        } 
      } 
    } 
  maf = aster[-1, ]
  maf = maf[order(maf$Batch,maf$Sample,maf$Chrom,maf$Start), ]

  # Save cleaned result.
  write.table(maf, paste0(lock.dir, '/', sam, '_clean.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
