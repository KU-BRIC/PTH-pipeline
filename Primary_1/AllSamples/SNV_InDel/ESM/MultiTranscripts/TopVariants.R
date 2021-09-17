####################################################################################################################
####################################################################################################################
# Create tables with variants with top ESM scores.
# Author: Haiying Kong
# Last Modified: 6 September 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

# Set parameters.
n.var = 500

####################################################################################################################
####################################################################################################################
# Read in all SAAS variants.
maf = read.table('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Apple_SAAS_full.txt',
                 header=TRUE, quote='', fill=TRUE, sep='\t')
N = length(unique(maf$Sample))

# Get ESM scores.
names(maf)[79:83] = paste0('ESM_', 1:5)
maf[ ,79:83] = apply(maf[ ,79:83], 2, function(x)  {
                                        sapply(x, function(y) max(as.numeric(unlist(strsplit(y, ',')))))
                                        })
maf$ESM_max = apply(maf[ ,79:83], 1, max)

####################################################################################################################
####################################################################################################################
# By sample_variant:
for (colname in paste0('ESM_', c(1:5, 'max')))  {
  if (nrow(maf)<=n.var)  {
    aster = maf[order(-abs(maf[ ,colname])), ]
    }    else  {
    aster = maf[order(-abs(maf[ ,colname]))[1:n.var], ]
    }
  write.table(aster,
              paste0('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/TopVariants_SAAS_', colname, '_bySampleVariant.maf'),
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
# By variant:
maf$N_Tumor = 1
formul = as.formula(paste0('N_Tumor~', paste(names(maf)[3:10], collapse='+')))
counts = aggregate(formul, data=maf, sum)
counts$N_Tumor_pct = round(counts$N_Tumor/N, 4)

idx = which(duplicated(maf[ ,3:10]))
maf = maf[-idx, -c(1:2,ncol(maf))]

for (colname in paste0('ESM_', c(1:5, 'max')))  {
  if (nrow(maf)<=n.var)  {
    aster = maf[order(-abs(maf[ ,colname])), ]
    }    else  {
    aster = maf[order(-abs(maf[ ,colname]))[1:n.var], ]
    }
  aster = merge(counts, aster, by=names(maf)[1:8])
  aster = aster[ ,c(1:8,match(colname,names(aster)),9:78)]
  aster = aster[order(-abs(aster[ ,colname])), ]
  write.table(aster, 
              paste0('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/TopVariants_SAAS_', colname, '_byVariant.maf'),
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
