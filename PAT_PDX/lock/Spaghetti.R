####################################################################################################################
####################################################################################################################
# Plot Spaghetti plots.
# Author: Haiying Kong
# Last Modified: 15 December 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(reshape)
library(ggplot2)

thresh = 0.05

# Clean the result directory.
res.dir = '/home/projects/cu_10184/projects/PTH/PAT_PDX/Result/Evolution/Spaghetti'
if (dir.exists(res.dir))  unlink(res.dir, recursive=TRUE)
dir.create(res.dir)

####################################################################################################################
####################################################################################################################
ped = read.table('/home/projects/cu_10184/projects/PTH/PAT_PDX/Meta/Pedigree.txt',
                 header=TRUE, sep='\t')

maf.cols = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'PAT_AF', 'PDX0_AF')

for (i in 1:nrow(ped))   {
  file.name = paste(ped$PAT_Batch[i], ped$PAT_Sample[i], ped$PDX0_Batch[i], ped$PDX0_Sample[i], sep='_')
  aster = read.table(paste0('/home/projects/cu_10184/projects/PTH/PAT_PDX/Result/Evolution/Table/', file.name, '.txt'),
                     header=TRUE, quote=' ', sep='\t')[ ,maf.cols]
  aster[ ,3:4] = apply(aster[ ,3:4], 2, as.character)
  aster$ID = apply(aster[ ,1:4], 1, function(x) paste(x, collapse='_'))
  aster = aster[ ,-(2:6)]

  # Filter out low AF variants.
  aster = aster[(aster$PAT_AF>thresh | aster$PDX0_AF>thresh), ]

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

  aster = aster[(aster$PAT_AF>thresh | aster$PDX0_AF>thresh), ]

  aster = melt(aster, id=c('Hugo_Symbol', 'ID'))
  names(aster)[3:4] = c('Generation', 'AF')
  aster$Generation = gsub('_AF', '', aster$Generation)
  aster$Gen = 1
  aster$Gen[aster$Generation=='PDX0'] = 2
  aster$Gen = as.integer(aster$Gen)

  # Plot:
  pdf(paste0(res.dir, '/', file.name, '.pdf'))
  p = ggplot(aster, aes(x=Gen, y=AF, color=factor(ID))) +
        geom_line() + geom_point() +
        scale_x_continuous(breaks=c(1,2), labels=c('1'='PAT', '2'='PDX0'), limits=c(1,2)) +
        scale_colour_discrete(name=' ') +
        theme(legend.text=element_text(size=1)) +
        theme_bw()
  print(p)
  dev.off()
  }

####################################################################################################################
####################################################################################################################
