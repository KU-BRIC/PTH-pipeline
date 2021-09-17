####################################################################################################################
####################################################################################################################
# Run ESM for one sample.
# Author: Haiying Kong
# Last Modified: 29 August 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

options(stringsAsFactors=FALSE)
rm(list=ls())

library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
hasProteinData(edb)

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
dir.name = args[1]
batch = args[2]
sam = args[3]

esm.dir = paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript')

####################################################################################################################
####################################################################################################################
# Define function to check single amino acid substitution.
amino.acids = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
check.saas = function(p.change)  {
               prt = sub('p.', '', p.change)
               ref = substr(prt, 1, 1)
               alt = substr(prt, nchar(prt), nchar(prt))
               loc = substr(prt, 2, nchar(prt)-1)
               flag = as.integer((ref %in% amino.acids) & (alt %in% amino.acids) & !is.na(suppressWarnings(as.integer(loc))))
               return(list(flag, prt, ref, alt, loc))
               }

trans = function(other.trans)  {
          flag = 0
          enst = c()
          prt = c()
          ref = c()
          alt = c()
          loc = c()
          x = unlist(strsplit(other.trans, '\\|'))
          idx = grep('_Missense_Mutation_p.', x)
          if (length(idx)>0) {
            x = x[idx]
            for (x1 in x)  {
              tmp = unlist(strsplit(x1, '_Missense_Mutation_p.'))
              tmp.prt = unlist(check.saas(tmp[2]))
              if (tmp.prt[1]==1)   {
                flag = 1
                enst = c(enst, paste0('ENST', unlist(strsplit(unlist(strsplit(tmp[1], '_ENST'))[2], '\\.'))[1]))
                prt = c(prt, tmp.prt[2])
                ref = c(ref, tmp.prt[3])
                alt = c(alt, tmp.prt[4])
                loc = c(loc, tmp.prt[5])
                }
              }
            }
          return(list(flag, enst, prt, ref, alt, loc))
          }

####################################################################################################################
####################################################################################################################
# Read in maf file.
maf = read.table(paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/TechError/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')

# idx for SAAS variants.
idx = which(maf$Protein_Change!='' & maf$Reference_Allele!='-' & maf$Tumor_Seq_Allele2!='-')

# Save non-SAAS.
write.table(maf[-idx, ], paste0(esm.dir, '/MultiTranscripts/', sam, '_rest.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# Save header for final SAAS file.
header = c(names(maf), paste0('esm1v_t33_650M_UR90S_', 1:5))
header = t(header)
write.table(header, paste0(esm.dir, '/MultiTranscripts/', sam, '_SAAS.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# Run ESM for candidates for SAAS.
maf = maf[idx, ]
for (i in 1:nrow(maf))   {
  flag = 0
  enst = c()
  prt = c()
  ref = c()
  alt = c()
  loc = c()

  # Protein_Change:
  tmp.prt = unlist(check.saas(maf$Protein_Change[i]))
  if (tmp.prt[1]==1)  {
    flag = 1
    enst = c(enst, unlist(strsplit(maf$Annotation_Transcript[i], '\\.'))[1])
    prt = c(prt, tmp.prt[2])
    ref = c(ref, tmp.prt[3])
    alt = c(alt, tmp.prt[4])
    loc = c(loc, tmp.prt[5])
    }

  # Other_Transcripts:
  other.trans = trans(maf$Other_Transcripts[i])


####################################################################################################################
####################################################################################################################


  prts = proteins(edb, filter=GeneNameFilter(maf$Hugo_Symbol[i]), return.type = "AAStringSet")
  idx.prt = which(mcols(prts)$tx_id == Annotation_Transcript[i])

  if (length(idx.prt)==0)  {
    write.table(maf[i, ], paste0(esm.dir, '/MultiTranscripts/', sam, '_rest.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
    next
    }

  prt.seq = as.data.frame(prts)[idx.prt, 1]
  prt.len = nchar(prt.seq)
  if (prt.len <= 1024)  {
    seq = prt.seq
    left = 1
    }  else  {
    left = max(1, aa.loc[i]-500)
    right = min(prt.len, aa.loc[i]+500)
    if ((right-left<1000) & (prt.len>1001))  {
      if (left==1)  right=min(prt.len, 1001)
      if (right==prt.len)  left=max(1, prt.len-1000)
      }
    seq = substr(prt.seq, left, right)
    }

  write.table(seq, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_seq.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  write.table(left, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_offset.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

  # maf file:
  write.table(maf[i, ], paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_acorn.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

  # input file.
  aster = as.data.frame(Protein_Change[i])
  names(aster) = 'Protein_Change'
  write.table(aster, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_input.csv'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=',')

  # Run ESM:
  system(paste0('sh /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/Filter/ESM/OneTranscript/esm_one_variant.sh -d ', dir.name, ' -b ', batch, ' -s ', sam, ' -i ', as.character(i)))
  }

####################################################################################################################
####################################################################################################################
