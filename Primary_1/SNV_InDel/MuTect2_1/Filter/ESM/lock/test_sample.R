####################################################################################################################
####################################################################################################################
# Run ESM for one sample.
# Author: Haiying Kong
# Last Modified: 25 August 2021
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

esm.dir = paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/ESM')

####################################################################################################################
####################################################################################################################
# Read in maf file.
maf = read.table(paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/TechError/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')
idx = which(maf$Variant_Type!='SNP' | maf$Variant_Classification!='Missense_Mutation' | maf$Protein_Change=='' | maf$Annotation_Transcript=='')
maf$Annotation_Transcript = sapply(maf$Annotation_Transcript, function(x)  unlist(strsplit(x, '\\.'))[1])

write.table(maf[idx, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

maf = maf[-idx, ]
Protein_Change = gsub('p.', '', maf$Protein_Change)

# Save header for final result file.
header = c(names(maf), paste0('esm1v_t33_650M_UR90S_', 1:5))
header = t(header)
write.table(header, paste0(esm.dir, '/Original/', sam, '_Missense.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

# Get amino acid locations.
aa.loc = sapply(maf$Protein_Change, function(x)  as.integer(substr(x,4,nchar(x)-1)))

for (i in 1:20)   {
  # protein sequence:
  prts = proteins(edb, filter=GeneNameFilter(maf$Hugo_Symbol[i]), return.type = "AAStringSet")
  idx.prt = which(mcols(prts)$tx_id == maf$Annotation_Transcript[i])

  if (length(idx.prt)==0)  {
    write.table(maf[i, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
    }  else  {
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
    system(paste0('sh /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/Filter/ESM/test_variant.sh -d ', dir.name, ' -b ', batch, ' -s ', sam, ' -i ', as.character(i)))
    }
  }

####################################################################################################################
####################################################################################################################
