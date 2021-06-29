####################################################################################################################
####################################################################################################################
# Get list of all variants from all callers and plot venn diagram.
# Author: Haiying Kong
# Last Modified: 7 June 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)
library(limma)
library(tidyverse)
library(ggforce)
library(gridExtra)

####################################################################################################################
####################################################################################################################
# Get batch list.
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

####################################################################################################################
# Before filtering:
####################################################################################################################
# Collect variants from all samples.
aster = c()
for (batch in batches)  {
  maf.files = dir(paste0('BatchWork/', batch, '/Result/SNV_InDel/AllVariants/Callers_Wide'), pattern='.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0('BatchWork/', batch, '/Result/SNV_InDel/AllVariants/Callers_Wide/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
    maf = maf[ ,c(1:6,8,10,13,12,11,24,19)]
    aster = rbind(aster, maf)
    }
  }
names(aster) = c('Batch', 'Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'AF_VarDict', 'AF_SNVer', 'AF_LoFreq', 'N_Alt', 'DP')

# Count number of variants filtered by sequence quality for each caller.
thresh.n.alt = 4
thresh.dp = 200
thresh.high.dp = read.table('/home/projects/cu_10184/projects/PTH/Reference/Filtering/Thresh_HighDP.txt',
                            header=FALSE, sep='\t')[1,1]
aster$Flag = 0
aster$Flag[(aster$N_Alt>=thresh.n.alt & aster$DP>=thresh.dp & aster$DP<=thresh.high.dp)] = 1
counts = data.frame(VarDict = c(sum(aster$AF_VarDict!=0), sum(aster$AF_VarDict!=0 & aster$Flag==0)),
                    SNVer = c(sum(aster$AF_SNVer!=0), sum(aster$AF_SNVer!=0 & aster$Flag==0)),
                    LoFreq = c(sum(aster$AF_LoFreq!=0), sum(aster$AF_LoFreq!=0 & aster$Flag==0)))
row.names(counts) = c('N_All', 'N_Filtered')
write.table(counts, 'AllBatches/SNV_InDel/CountVariants_3Callers_Before_After_Filter.txt', row.names=TRUE, col.names=TRUE, quote=FALSE, sep='\t')

aster = aster[ ,1:11]

# Plot venn diagram.
aster$Tag = apply(aster[ ,1:8], 1, function(x) paste(x, collapse='_'))
aster = aster[ ,c(12,9:11)]
apple = data.frame(VarDict = (aster$AF_VarDict!=0),
                   SNVer = (aster$AF_SNVer!=0),
                   LoFreq = (aster$AF_LoFreq!=0))

df.venn = data.frame(x = c(0, 0.866, -0.866),
                     y = c(1, -0.5, -0.5),
                     labels = c('VarDict', 'SNVer', 'LoFreq'))
vdc = vennCounts(apple)
class(vdc) = 'matrix'
df.vdc = as.data.frame(vdc)[-1,] %>% mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0), y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))

all = ggplot(df.venn) +
        ggtitle('Before Filtering') +
        geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .3, size = 1, colour = 'grey') +
        coord_fixed() +
        theme_void() +
        theme(legend.position = 'bottom') +
        scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
        scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
        labs(fill = NULL) +
        annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5)

####################################################################################################################
# After filtering.
####################################################################################################################
for (filter.scheme in c('Long', 'Medium', 'Short'))  {

  # Collect variants from all samples.
  aster = c()
  for (batch in batches)  {
    maf.files = dir(paste0('BatchWork/', batch, '/Result/SNV_InDel/Filtered/', filter.scheme), pattern='.maf')
    for (maf.file in maf.files)  {
      maf = read.table(paste0('BatchWork/', batch, '/Result/SNV_InDel/Filtered/', filter.scheme, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
      maf = maf[ ,c(1:6,8,10,13,12,11)]
      aster = rbind(aster, maf)
      }
    }
  names(aster) = c('Batch', 'Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'AF_VarDict', 'AF_SNVer', 'AF_LoFreq')

  # Plot venn diagram.
  aster$Tag = apply(aster[ ,1:8], 1, function(x) paste(x, collapse='_'))
  aster = aster[ ,c(12,9:11)]
  apple = data.frame(VarDict = (aster$AF_VarDict!=0),
                     SNVer = (aster$AF_SNVer!=0),
                     LoFreq = (aster$AF_LoFreq!=0))
                   

  df.venn = data.frame(x = c(0, 0.866, -0.866),
                       y = c(1, -0.5, -0.5),
                       labels = c('VarDict', 'SNVer', 'LoFreq'))
  vdc = vennCounts(apple)
  class(vdc) = 'matrix'
  df.vdc = as.data.frame(vdc)[-1,] %>% mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0), y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))

  assign(filter.scheme,
    ggplot(df.venn) +
      ggtitle(paste0('Filtered by ', filter.scheme, ' Scheme')) +
      geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .3, size = 1, colour = 'grey') +
      coord_fixed() +
      theme_void() +
      theme(legend.position = 'bottom') +
      scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
      scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
      labs(fill = NULL) +
      annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5))
  }

####################################################################################################################
pdf('AllBatches/SNV_InDel/vennDiagram_3Callers_Before_After_Filtering.pdf')

margin = theme(plot.margin = unit(rep(0.5, 4), "cm"))
grid.arrange(all, Long, Medium, Short, nrow=2, ncol=2)

dev.off()

####################################################################################################################
####################################################################################################################
