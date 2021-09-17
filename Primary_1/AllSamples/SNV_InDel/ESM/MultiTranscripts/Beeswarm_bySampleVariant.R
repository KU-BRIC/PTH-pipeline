####################################################################################################################
####################################################################################################################
# 
# Author: Haiying Kong
# Last Modified: 5 September 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(data.table)
library(beeswarm)

####################################################################################################################
####################################################################################################################
#
####################################################################################################################
# Read in table with all SAAS and trim it.
maf = read.table('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Apple_SAAS.txt',
                 header=TRUE, quote='', fill=TRUE, sep='\t')
names(maf)[20:24] = paste0('ESM_', 1:5)
maf[ ,20:24] = apply(maf[ ,20:24], 2, function(x)  {
                                        sapply(x, function(y) max(as.numeric(unlist(strsplit(y, ',')))))
                                        })
maf = as.data.frame(maf)
maf$ESM_max = apply(maf[ ,20:24], 1, max)
maf$Patient = as.vector(sapply(maf$Sample, function(x)  unlist(strsplit(x, '-'))[1]))
maf = maf[ ,c(1,26,2,3,13,16,17,20:25)]

####################################################################################################################
# With Horizen and REDCap:
####################################################################################################################
# Read in known Horizon and REDCap variants.
tag = read.table('/home/projects/cu_10184/projects/PTH/AllBatches/Result/SNV_InDel/ClinicTag/CalledReferenceVariants.txt',
                 header=TRUE, sep='\t')
tag = tag[ ,c('Batch', 'PTH_ID', 'Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'Ref', 'Alt', 'Refseq_ID', 'cDNA_Change', 'Protein_Change', 'AF_Clinic', 'Note')]

Protein_Change = gsub('p.', '', tag$Protein_Change)
amino.acids = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
flag = sapply(Protein_Change, function(x)  {
                                ref = substr(x, 1, 1)
                                alt = substr(x, nchar(x), nchar(x))
                                loc = substr(x, 2, nchar(x)-1)
                                flag = as.integer((ref %in% amino.acids) & (alt %in% amino.acids) & !is.na(suppressWarnings(as.integer(loc))))
                                })
idx = which(as.vector(flag==1))
tag = tag[idx, ]
tag.label = paste(tag$PTH_ID, tag$Hugo_Symbol, tag$Protein_Change, sep='_')

# Tag our variant list.
aster = maf[ ,c(1:5,8:13)]
aster.label = paste(aster$Patient, aster$Hugo_Symbol, aster$Protein_Change, sep='_')

idx = match(tag.label, aster.label)
idx = idx[!is.na(idx)]

aster$color = 'green'
aster$color[idx] = 'red'

ast = aster[ ,c(12,6:11)]
ast = melt(setDT(ast), id.vars='color', variable.name='ESM_ID')
names(ast)[3] = 'ESM_Score'

####################################################################################################################
# Plot Beeswarm.
pdf('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Beeswarm_ESM_HorizenREDCap_bySampleVariant.pdf')

# All ESM scores.
beeswarm(ESM_Score ~ ESM_ID, data=ast, method='swarm', pch=16, cex=0.15, pwcol=color,
         xlab='ESM_ID', ylab='ESM_Score', labels = c(paste0('ESM_',1:5),'ESM_max'))
legend('topright', legend=c('Horizon_REDCap','Not'), pch=16, col=c('red','green'))

# By each ESM scores.
for (j in 6:11)  {
  beeswarm(aster[ ,j], method='swarm', pch=16, cex=0.2, pwcol=aster$color,
           xlab=names(aster)[j], ylab='ESM_Score')
  legend('topright', legend=c('Horizon_REDCap','Not'), pch=16, col=c('red','green'))
  }

dev.off()

####################################################################################################################
# With COSMIC_tissue_types_affected:
####################################################################################################################
# Get color for variants.
aster = maf[ ,-(6:7)]

aster$color = 'green'
aster$color[maf$COSMIC_tissue_types_affected!=''] = 'red'

ast = aster[ ,c(12,6:11)]
ast = melt(setDT(ast), id.vars='color', variable.name='ESM_ID')
names(ast)[3] = 'ESM_Score'

####################################################################################################################
# Plot Beeswarm.
pdf('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Beeswarm_ESM_COSMIC_long_bySampleVariant.pdf')

# All ESM scores.
beeswarm(ESM_Score ~ ESM_ID, data=ast, method='swarm', pch=16, cex=0.15, pwcol=color,
         xlab='ESM_ID', ylab='ESM_Score', labels = c(paste0('ESM_',1:5),'ESM_max'))
legend('topright', legend=c('COSMIC_long','Not'), pch=16, col=c('red','green'))

# By each ESM scores.
for (j in 6:11)  {
  beeswarm(aster[ ,j], method='swarm', pch=16, cex=0.2, pwcol=aster$color,
           xlab=names(aster)[j], ylab='ESM_Score')
  legend('topright', legend=c('COSMIC_long','Not'), pch=16, col=c('red','green'))
  }

dev.off()

####################################################################################################################
# With COSMIC_overlapping_mutations:
####################################################################################################################
# Get color for variants.
aster = maf[ ,-(6:7)]

aster$color = 'green'
aster$color[maf$COSMIC_overlapping_mutations!=''] = 'red'

ast = aster[ ,c(12,6:11)]
ast = melt(setDT(ast), id.vars='color', variable.name='ESM_ID')
names(ast)[3] = 'ESM_Score'

####################################################################################################################
# Plot Beeswarm.
pdf('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Beeswarm_ESM_COSMIC_short_bySampleVariant.pdf')

# All ESM scores.
beeswarm(ESM_Score ~ ESM_ID, data=ast, method='swarm', pch=16, cex=0.15, pwcol=color,
         xlab='ESM_ID', ylab='ESM_Score', labels = c(paste0('ESM_',1:5),'ESM_max'))
legend('topright', legend=c('COSMIC_short','Not'), pch=16, col=c('red','green'))

# By each ESM scores.
for (j in 6:11)  {
  beeswarm(aster[ ,j], method='swarm', pch=16, cex=0.2, pwcol=aster$color,
           xlab=names(aster)[j], ylab='ESM_Score')
  legend('topright', legend=c('COSMIC_short','Not'), pch=16, col=c('red','green'))
  }

dev.off()

####################################################################################################################
####################################################################################################################
