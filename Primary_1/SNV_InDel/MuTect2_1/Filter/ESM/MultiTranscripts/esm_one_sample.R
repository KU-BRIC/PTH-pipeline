####################################################################################################################
####################################################################################################################
# Run ESM for one sample.
# Author: Haiying Kong
# Last Modified: 30 August 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

options(stringsAsFactors=FALSE)
rm(list=ls())

library(biomaRt)
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")


####################################################################################################################
# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
dir.name = args[1]
batch = args[2]
sam = args[3]

esm.dir = paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts')

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
write.table(maf[-idx, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# Exclude non-SAAS.
maf = maf[idx, ]

# Save header for final SAAS file.
header = c(names(maf), paste0('esm1v_t33_650M_UR90S_', 1:5))
header = t(header)
write.table(header, paste0(esm.dir, '/Original/', sam, '_SAAS.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# Process SAAS candidates.
for (i in 1:nrow(maf))  {

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
  if (other.trans[[1]]==1)  {
    enst = c(enst, other.trans[[2]])
    prt = c(prt, other.trans[[3]])
    ref = c(ref, other.trans[[4]])
    alt = c(alt, other.trans[[5]])
    loc = c(loc, other.trans[[6]])
    }

  if (flag==0)  {
    write.table(maf[i, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
    next
    }

  # Get peptide sequence for the transcripts.
  seqs = getSequence(id=enst, seqType="peptide", mart=mart, type="ensembl_transcript_id")

  if (nrow(seqs)==0)  {
    write.table(maf[i, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
    next
    }

  idx = match(seqs$ensembl_transcript_id, enst)

  # Get all SAAS candidates for the variant.
  aster = data.frame(Seq = seqs$peptide,
                     Protein_Change = prt[idx],
                     Ref = ref[idx],
                     Loc = as.numeric(loc)[idx]
                    )
  aster = unique(aster)

  if (nrow(aster)==0)  {
    write.table(maf[i, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
    next
    }

  write.table(c(), paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_esm.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

  i.flag = 0
  for (j in 1:nrow(aster))  {
    if (substr(aster$Seq[j], aster$Loc[j], aster$Loc[j]) != aster$Ref[j])  next  
    i.flag = 1
    print(i.flag)
    prt.len = nchar(aster$Seq[j])
    if (prt.len <= 1024)  {
      seq = aster$Seq[j]
      left = 1
      }  else  {
      left = max(1, aster$Loc[j]-500)
      right = min(prt.len, aster$Loc[j]+500)
      if ((right-left<1000) & (prt.len>1001))  {
        if (left==1)  right=min(prt.len, 1001) 
        if (right==prt.len)  left=max(1, prt.len-1000)
        }
      seq = substr(aster$Seq[j], left, right)
      }
    
    write.table(seq, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_', as.character(j), '_seq.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    write.table(left, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_', as.character(j), '_offset.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  
    # input file.
    input = as.data.frame(aster$Protein_Change[j])
    names(input) = 'Protein_Change'
    write.table(input, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_', as.character(j), '_input.csv'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=',')
  
    # Run ESM:
    system(paste0('sh /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/Filter/ESM/MultiTranscripts/esm_one_seq.sh -d ', dir.name, ' -b ', batch, ' -s ', sam, ' -i ', as.character(i), ' -j ', as.character(j)))
    }

  if (i.flag==0)
    write.table(maf[i, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')    else  {
    one.var = read.table(paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_esm.txt'), header=FALSE, sep=',')
    one.var = apply(one.var, 2, function(x) paste(x, collapse=','))
    one.var = cbind(maf[i, ], t(one.var))
    write.table(one.var, paste0(esm.dir, '/Original/', sam, '_SAAS.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
    }

  system(paste0('rm -rf ', esm.dir, '/temp/', sam, '_', as.character(i), '_esm.txt'))

  }

####################################################################################################################
####################################################################################################################
