####################################################################################################################
####################################################################################################################
# 
# Author: Haiying Kong
# Last Modified: 15 December 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(parallel)
library(limma)
library(ggplot2)
library(ggforce)
library(gridExtra)
library(reshape)

# Set parameters.
filter.scheme = 'Long'
maf.cols = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'AF', 'DP')

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
thresh.spa = args[1]

# Clean lock and result directories.
lock.dir = '/home/projects/cu_10184/projects/PTH/PAT_PDX/Lock/SNV_InDel/maf'
if (dir.exists(lock.dir))  unlink(lock.dir, recursive=TRUE)
dir.create(lock.dir)
dir.create(paste0(lock.dir, '/PAT'))

res.dir = '/home/projects/cu_10184/projects/PTH/PAT_PDX/Result/Evolution'
if (dir.exists(res.dir))  unlink(res.dir, recursive=TRUE)
dir.create(res.dir)
dir.create(paste0(res.dir, '/Table'))
dir.create(paste0(res.dir, '/Plot'))

####################################################################################################################
# Read in patient samples, xenograft, and small panel target information.
pat.info = read.table('/home/projects/cu_10184/projects/PTH/Meta/SampleInfo.txt', header=TRUE, sep='\t')
pdx.info = read.table('/home/projects/cu_10184/projects/PTH/Meta/PDX_Info.txt', header=TRUE, sep='\t')

pat.pdx.info = read.table('/home/projects/cu_10184/projects/PTH/Meta/PAT_PDX_Info.txt', header=TRUE, sep='\t')
tmp = pat.pdx.info$Subportion_ID
tmp = sapply(tmp, function(x) unlist(strsplit(x, ' '))[1])
names(tmp) = NULL
tmp = sapply(tmp, function(x) unlist(strsplit(x, '-'))[1])
names(tmp) = NULL
pat.pdx.info$Subport_ID = tmp

target = read.table('/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/Focused_myeloid_panel-All_target_segments_covered_by_probes-TE-93310852_hg38_v2_190722165759.bed',
                    header=FALSE, sep='\t')
names(target) = c('Chrom', 'Start', 'End')

####################################################################################################################
# Trim pat.pdx.info file for samples we have.
sam.ids = unique(c(pat.info$PTH_ID, pdx.info$Mouse_ID))
pat.pdx.info = pat.pdx.info[(pat.pdx.info$Mouse_ID %in% sam.ids) & (pat.pdx.info$Subport_ID %in% sam.ids), ]

####################################################################################################################
# Create pedigree table.
ped = c()
pats = sort(unique(pat.pdx.info$Subport_ID[pat.pdx.info$Mouse_Passage=='P0']))

for (pat in pats)  {
  pat.idxs = which(pat.info$PTH_ID==pat)
  pdx0s = pat.pdx.info$Mouse_ID[pat.pdx.info$Subport_ID==pat & pat.pdx.info$Mouse_Passage=='P0']
  if (length(pat.idxs)==0 | length(pdx0s)==0) next

  for (pat.idx in pat.idxs)  {
    # Patient sample information:
    pat.batch = pat.info$Batch_ID[pat.idx]
    pat.pat = pat.info$PTH_ID[pat.idx]
    pat.sam = pat.info$Sample_ID[pat.idx]

    # PDX0 information:
    for (pdx0 in pdx0s)  {
      pdx0.idxs = which(pdx.info$Mouse_ID==pdx0)
      for (pdx0.idx in pdx0.idxs)  {
        pdx0.batch = pdx.info$Batch_ID[pdx0.idx]
        pdx0.mouse = pdx.info$Mouse_ID[pdx0.idx]
        pdx0.sam = pdx.info$Sample_ID[pdx0.idx]

        aster = c(pat.batch, pat.pat, pat.sam, pdx0.batch, pdx0.mouse, pdx0.sam, '', '', '')

        # Check for PDX1.
        pdx1s = pat.pdx.info$Mouse_ID[pat.pdx.info$Subport_ID==pdx0.mouse & pat.pdx.info$Mouse_Passage=='P1']
        if (length(pdx1s)>0)  {
          for (pdx1 in pdx1s)  {
            pdx1.idxs = which(pdx.info$Mouse_ID==pdx1)
            for (pdx1.idx in pdx1.idxs)  {
              pdx1.batch = pdx.info$Batch_ID[pdx1.idx]
              pdx1.mouse = pdx.info$Mouse_ID[pdx1.idx]
              pdx1.sam = pdx.info$Sample_ID[pdx1.idx]
              aster = c(pat.batch, pat.pat, pat.sam, pdx0.batch, pdx0.mouse, pdx0.sam, pdx1.batch, pdx1.mouse, pdx1.sam)
              }
            }
          }
        ped = rbind(ped, aster)
        }
      }
    }
  }

ped = as.data.frame(ped)
row.names(ped) = NULL
names(ped) = c('PAT_Batch', 'PAT_Patient', 'PAT_Sample', 'PDX0_Batch', 'PDX0_Mouse', 'PDX0_Sample', 'PDX1_Batch', 'PDX1_Mouse', 'PDX1_Sample')

write.table(ped, '/home/projects/cu_10184/projects/PTH/PAT_PDX/Meta/Pedigree.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
# venn Diagram count and plot, spaghetti plot.
####################################################################################################################

venn = c()

for (i in 1:nrow(ped))   {
  ##############################################################
  # Variants in PAT:
  # vcf:
  pat.vcf = read.table(paste0('/home/projects/cu_10184/projects/PTH/BatchWork/', ped$PAT_Batch[i], '/Lock/SNV_InDel/VarDict/vcf/', ped$PAT_Sample[i], '.vcf'),
                       header=FALSE, quote='', sep='\t')

  # Filter for the smaller target regions.
  n.cores = detectCores()
  ans = mclapply(1:nrow(pat.vcf), function(ii)  {
                                    idx = which(target$Chrom==pat.vcf$V1[ii] & target$Start<=pat.vcf$V2[ii] & target$End>=(pat.vcf$V2[ii]+max(nchar(pat.vcf$V4[ii]),nchar(pat.vcf$V5[ii]))))
                                    if (length(idx)>0)
                                      flag = 1   else
                                      flag = 0
                                    flag
                                    },
                 mc.cores = n.cores)
  pat.vcf = pat.vcf[ans==1, ]
  write.table(pat.vcf, paste0('/home/projects/cu_10184/projects/PTH/PAT_PDX/Lock/SNV_InDel/vcf/PAT/', ped$PAT_Sample[i], '.vcf'),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

  # maf:
  pat.maf = read.table(paste0('/home/projects/cu_10184/projects/PTH/BatchWork/', ped$PAT_Batch[i], '/Result/SNV_InDel/Filtered/', filter.scheme, '/', ped$PAT_Sample[i], '.maf'),
                       header=TRUE, quote='', sep='\t')

  # Filter for the smaller target regions.
  n.cores = detectCores()
  ans = mclapply(1:nrow(pat.maf), function(i)  {
                                    idx = which(target$Chrom==pat.maf$Chromosome[i] & target$Start<=pat.maf$Start_Position[i] & target$End>=pat.maf$End_Position[i])
                                    if (length(idx)>0)
                                      flag = 1   else
                                      flag = 0
                                    flag
                                    },
                 mc.cores = n.cores)
  pat.maf = pat.maf[ans==1, ]
  write.table(pat.maf, paste0('/home/projects/cu_10184/projects/PTH/PAT_PDX/Lock/SNV_InDel/vcf/PAT/', ped$PAT_Sample[i], '.maf'),
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  pat = pat.maf

  ##############################################################
  # Variants in PDX0:
  pdx0 = read.table(paste0('/home/projects/cu_10184/projects/PTH/BatchWork/', ped$PDX0_Batch[i], '/Result/SNV_InDel/Filtered/', filter.scheme, '/', ped$PDX0_Sample[i], '.maf'),
                    header=TRUE, quote='', sep='\t')[ ,1:117]

  ##############################################################
  # Combine variants in PAT and PDX0:
  apple = c()

  if (nrow(pat)>0)  pat[ ,5:6] = apply(pat[ ,5:6], 2, as.character)
  pat.flag = apply(pat[ ,maf.cols[1:6]], 1, function(x) paste(x, collapse='_'))
  if (nrow(pat)>0)  pat[ ,5:6] = apply(pat[ ,5:6], 2, as.numeric)

  if (nrow(pdx0)>0)  pdx0[ ,5:6] = apply(pdx0[ ,5:6], 2, as.character)
  pdx0.flag = apply(pdx0[ ,maf.cols[1:6]], 1, function(x) paste(x, collapse='_'))
  if (nrow(pdx0)>0)  pdx0[ ,5:6] = apply(pdx0[ ,5:6], 2, as.numeric)

  idx = match(pat.flag, pdx0.flag)

  # Variants in both pat and pdx0.
  if (length(which(!is.na(idx))) > 0)  {
    apple = cbind(pat[(!is.na(idx)), ], pdx0[idx[!is.na(idx)], c('Batch', 'Sample', 'AF', 'DP')])
    apple = apple[ ,c(1:2,118:119,3:19,120:121,20:117)]
    names(apple)[c(1:4, 20:23)] = c('PAT_Batch', 'PAT_Sample', 'PDX0_Batch', 'PDX0_Sample', 'PAT_AF', 'PAT_DP', 'PDX0_AF', 'PDX0_DP')
    }

  # Variants in only pat.
  if (length(which(is.na(idx))) > 0)  {
    aster = pat[is.na(idx), ]
    aster = cbind(aster[ ,1:2], pdx0$Batch[1], pdx0$Sample[1], aster[ ,3:19], 0.0, 0, aster[ ,20:117])
    names(aster)[c(1:4, 20:23)] = c('PAT_Batch', 'PAT_Sample', 'PDX0_Batch', 'PDX0_Sample', 'PAT_AF', 'PAT_DP', 'PDX0_AF', 'PDX0_DP')
    if (is.null(apple))
      apple = aster  else
      apple = rbind(apple, aster)
    }

  # Variants in only pdx0.
  idx.pdx0 = which(!(pdx0.flag %in% pat.flag))
  if (length(idx.pdx0) > 0)  {
    aster = pdx0[idx.pdx0, ]
    aster = cbind(pat$Batch[1], pat$Sample[1], aster[ ,1:17], 0.0, 0, aster[ ,18:117])
    names(aster)[c(1:4, 20:23)] = c('PAT_Batch', 'PAT_Sample', 'PDX0_Batch', 'PDX0_Sample', 'PAT_AF', 'PAT_DP', 'PDX0_AF', 'PDX0_DP')
    if (is.null(apple))
      apple = aster  else
      apple = rbind(apple, aster)
    }

  apple$Delta_AF = apple$PDX0_AF - apple$PAT_AF
  apple = apple[ ,c(1:23,122,24:121)]

  ##############################################################
  # Save apple.
  write.table(apple, paste0(res.dir, '/Table/', ped$PAT_Batch[i], '_', ped$PAT_Sample[i], '_', ped$PDX0_Batch[i], '_', ped$PDX0_Sample[i], '.txt'),
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

  ##############################################################
  # Count variants in PAT and PDX0, and update venn.
  ven = c(length(which(apple$PAT_DP>0)), length(which(apple$PDX0_DP>0)), length(which(apple$PAT_DP>0 & apple$PDX0_DP>0)))
  venn = rbind(venn, c(ped$PAT_Batch[i], ped$PAT_Sample[i], ped$PDX0_Batch[i], ped$PDX0_Sample[i], ven))

  ##############################################################
  # Plot venn Diagram and spaghetti.
  ##############################################################
  file.name = paste(ped$PAT_Batch[i], ped$PAT_Sample[i], ped$PDX0_Batch[i], ped$PDX0_Sample[i], sep='_')
  pdf(paste0(res.dir, '/Plot/', file.name, '.pdf'))

  ###############
  # venn Diagram:
  aster = data.frame(PAT = (apple$PAT_DP!=0),
                     PDX0 = (apple$PDX0_DP!=0))
  vdc = vennCounts(aster)
  vennDiagram(vdc, main=paste0('\n', '\n', '\n', '\n', ped$PAT_Sample[i], '\n', ped$PDX0_Sample[i]), circle.col=c('cornflowerblue', 'firebrick'),
              cex.main=0.7, cex=0.7)

  ###############
  # Spaghetti:
  maf.cols = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Protein_Change', 'PAT_AF', 'PDX0_AF')
  aster = apple[ ,maf.cols]

  # Filter out low AF variants.
  aster = aster[(aster$PAT_AF>thresh.spa | aster$PDX0_AF>thresh.spa), ]

  aster[ ,3:4] = apply(aster[ ,3:4], 2, as.character)
  for (i in 1:nrow(aster))  {
    if (nchar(aster$Protein_Change[i])>12)
    aster$Protein_Change[i] = substr(aster$Protein_Change[i], 1, 12)
    }

  aster$ID = apply(aster[ ,c(1:4,7)], 1, function(x) paste(x, collapse='_'))
  aster = aster[ ,-(1:7)] 
  
  # Treat duplicated IDs.
  idx = which(duplicated(aster$ID))
  idx.i = 0
  while (length(idx) >0)  {
    idx.i = idx.i + 1
    for (idx.idx in idx)  {
      aster$ID[idx.idx] = paste(unlist(strsplit(aster$ID[idx.idx], '\\.'))[1], idx.i, sep='\\.')
      }
    idx = which(duplicated(aster$ID))
    }
    
  # Reshape:
  aster = melt(aster, id=c('ID'))
  names(aster)[2:3] = c('Generation', 'AF')
  aster$Generation = gsub('_AF', '', aster$Generation)
  aster$Gen = 1
  aster$Gen[aster$Generation=='PDX0'] = 2
  aster$Gen = as.integer(aster$Gen)
  
  # Plot:
  p = ggplot(aster, aes(x=Gen, y=AF, color=factor(ID))) +
        ggtitle(paste0(ped$PAT_Sample[i], '\n', ped$PDX0_Sample[i])) +
        labs(x='') +
        geom_line() + geom_point() + 
        scale_x_continuous(breaks=c(1,2), labels=c('1'='PAT', '2'='PDX0'), limits=c(1,2)) +
        scale_colour_discrete(name=' ') + 
        theme(plot.title=element_text(size=0.1, face='bold'), legend.text=element_text(size=0.7)) +
        theme_bw()
  print(p)

  dev.off()
  }

venn = as.data.frame(venn)
row.names(venn) = NULL
names(venn) = c('PAT_Batch', 'PAT_Sample', 'PDX0_Batch', 'PDX0_Sample', 'N_PAT', 'N_PDX0', 'N_Both')
venn[ ,5:7] = apply(venn[ ,5:7], 2, as.numeric)
write.table(venn, '/home/projects/cu_10184/projects/PTH/PAT_PDX/Result/Evolution/VariantCounts.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
