####################################################################################################################
####################################################################################################################
# Identify potential SNPs identified from PoN and region specific technical errors.
# Author: Haiying Kong
# Last Modified: 13 July 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10145/people/haikon/Project/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(parallel)
library(clValid)
library(reshape2)

####################################################################################################################
# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
scheme_dir = args[1]
scheme_name = args[2]

####################################################################################################################
# Define directory to save results.
ref.dir = paste0('/home/projects/cu_10184/projects/', scheme_dir, '/Reference/Filtering/', scheme_name)

# Read in thresholds.
thresh = read.table(paste0(ref.dir, '/Thresholds.txt'), header=TRUE, sep='\t')

####################################################################################################################
####################################################################################################################
# Read in sample information.
sam.info = read.table('Meta/SampleInfo.txt', header=TRUE, sep='\t')[ ,1:5]
sam.info = sam.info[order(sam.info$Batch_ID, sam.info$Sample_ID), ]

normals = unique(sam.info$Sample[sam.info$Group=='NORMAL'])
normals = sort(unique(sapply(normals, function(x) strsplit(x,'-')[[1]][1])))

####################################################################################################################
####################################################################################################################
# Get batch names.
batch.dir = paste0('/home/projects/cu_10184/projects/', scheme_dir, '/BatchWork/')
batches = dir(batch.dir, pattern='^Primary_')

####################################################################################################################
# Collect all DPs and AFs from all samples.
apple = c()
col.classes = c('NULL', rep('character',3), rep('integer',2), 'NULL', 'character', 'NULL', 'character', rep('NULL',7), 'numeric', 'integer', rep('NULL',98))

for (batch in batches)   {
  sam.dir = paste0(batch.dir, batch, '/Result/SNV_InDel/AllVariants/Callers_Wide/')
  samples = dir(sam.dir, pattern='.maf')
  for (sam in samples)   {
    aster = read.table(paste0(sam.dir, sam), header=TRUE, colClasses=col.classes, quote='', fill=TRUE, sep='\t')
    apple = rbind(apple, aster)
    }
  }
names(apple) = c('Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'AF', 'DP')

####################################################################################################################
# Create wide tables for AF, DP and Flag.
AF = apple[ ,c('Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'AF')]
formul = as.formula(paste0(paste(names(AF)[2:7], collapse='+'), ' ~ Sample'))
AF = dcast(AF, formul, value.var='AF')
AF[ ,-(1:6)][is.na(AF[ ,-(1:6)])] = 0

DP = apple[ ,c('Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'DP')]
formul = as.formula(paste0(paste(names(DP)[2:7], collapse='+'), ' ~ Sample'))
DP = dcast(DP, formul, value.var='DP')
DP[ ,-(1:6)][is.na(DP[ ,-(1:6)])] = 0

Flag = AF
Flag[ ,-(1:6)] = apply(Flag[ ,-(1:6)], 2, function(x)  {
                                            idx = which(x==0)
                                            res = rep(1, length(x))
                                            res[idx] = 0
                                            res   })

####################################################################################################################
# Find column index for normal and patients, and reorder the columns.
pth.id = sapply(names(AF)[-(1:6)], function(x) unlist(strsplit(x,'-'))[1])
names(pth.id) = NULL

idx.norm = which(pth.id %in% normals)
idx.norm = idx.norm[!is.na(idx.norm)]
idx.pat = (1:(ncol(AF)-6))[-idx.norm]
idx.norm = idx.norm + 6
idx.pat = idx.pat + 6

idx = c(1:6, idx.norm, idx.pat)

AF = AF[ ,idx]
DP = DP[ ,idx]
Flag = Flag[ ,idx]

pth.id = sapply(names(AF)[-(1:6)], function(x) unlist(strsplit(x,'-'))[1])
names(pth.id) = NULL

idx.norm = which(pth.id %in% normals)
idx.norm = idx.norm[!is.na(idx.norm)]
idx.pat = (1:(ncol(AF)-6))[-idx.norm]
idx.norm = idx.norm + 6
idx.pat = idx.pat + 6

####################################################################################################################
####################################################################################################################
# Filter with PoN.
####################################################################################################################
# Identify potential SNPs with PoN.
pth.id.norm = pth.id[idx.norm-6]
N = length(unique(pth.id.norm))
maf.pon = apply(Flag[ ,idx.norm], 1, function(x)   length(unique(pth.id.norm[which(x==1)]))) / N
idx = which(maf.pon > thresh$thresh_maf_norm)
SNPs = AF[idx, 1:6]

write.table(SNPs, paste0(ref.dir, '/SNPs_PoN.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# Filter with PoN.
AF = AF[-idx, ]
DP = DP[-idx, ]
Flag = Flag[-idx, ]

####################################################################################################################
####################################################################################################################
# Identify region specific technical errors with clustering.
####################################################################################################################
n.cores = detectCores()
ans = mclapply(1:nrow(AF), function(i)  {
                             af = as.numeric(t(AF[i,-(1:6)]))
                             flag = c(0, 0, 0)
                             
                             # Get information for this row.
                             idx.non0 = which(af>0)
                             af.non0 = as.data.frame(af[idx.non0])
                             row.names(af.non0) = names(AF[ ,-(1:6)])[idx.non0]
                             names(af.non0) = 'AF'
                             af.pat = AF[i,idx.pat]
                             af.norm = -(sort(-AF[i,idx.norm]))
                             
                             # Test for region specific technical error if 5 or more samples carry the variant.
                             if (nrow(af.non0)>=5 & sd(af.non0$AF)!=0)  {
                               # Check if AFs for this variants can be clearly clustered into 2.
                               score = clValid(af.non0, 2, clMethods="kmeans", validation="internal")
                               silhouette = as.numeric(measures(score,"Silhouette"))
                                   
                               if (silhouette>thresh$thresh_silhouette)  {
                                 # k-means clustering.
                                 km = kmeans(af.non0, 2)
                                 centers = km$centers 
                                 km.out = data.frame(Sample = row.names(af.non0),
                                                     AF = af.non0,
                                                     Cluster = km$cluster)
                                 row.names(km.out) = NULL
                                 if (centers[1] < centers[2])  {
                                   km.out$Cluster = -(km.out$Cluster - 2) + 1
                                   centers = c(centers[2], centers[1]) 
                                   }
                                    
                                 # Collect information for region specific technical error if the center of lower cluster is below threshold.
                                 if (centers[2]<thresh$thresh_lower_cluster_center)  {
                                   flag = c(1, centers[1], centers[2])
                                   }
                                 } 
                               } 
                             list(flag)
                             },
                mc.cores = n.cores)

flag = sapply(ans, function(x) list(x[[1]][[1]]))
centers.1 = round(unlist(sapply(ans, function(x) list(x[[1]][[2]]))), 4)
centers.2 = round(unlist(sapply(ans, function(x) list(x[[1]][[3]]))), 4)

idx = which(flag != 0)
apple = cbind(AF[idx, 1:6], centers.1[idx], centers.2[idx])
names(apple)[7:8] = c('UpperCenter', 'LowerCenter')

write.table(apple, paste0(ref.dir, '/RegionSpecificTechnicalError.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
