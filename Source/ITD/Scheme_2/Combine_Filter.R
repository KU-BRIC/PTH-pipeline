####################################################################################################################
####################################################################################################################
# Combine ITD results from VarDict, Pindel, ScanITD, getITD.
# Author: Haiying Kong
# Last Modified: 25 June 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set working directory, options and clean the space.
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
batch = args[1]
sam = args[2]
lock.dir = args[3]
res.dir = args[4]
thresh.n.alt = args[5]

####################################################################################################################
####################################################################################################################
# Read in ITD calls from tools and combine them.
####################################################################################################################

apple = c()

####################################################################################################################
# VarDict:
aster.file = paste0(lock.dir, '/VarDict/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  aster = aster[ ,c(names(aster)[1:9],'TYPE','cDNA_Change','Refseq_mRNA_Id','Codon_Change','Protein_Change')]
  names(aster)[10] = 'Type'
  aster = aster[aster$N_Alt>thresh.n.alt, ]
  if (nrow(aster)>0)  {
    aster = cbind('VarDict', aster)
    names(aster)[1] = 'Caller'
    apple = rbind(apple, aster)
    }
  }

####################################################################################################################
# Pindel:
aster.file = paste0(lock.dir, '/Pindel/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  if (nrow(aster)>0)  {
    aster = aster[aster$N_Alt>thresh.n.alt, ]
    if (nrow(aster)>0)  {
      aster = cbind('Pindel', aster)
      names(aster)[1] = 'Caller'
      aster$Type = paste0(aster$Type, ':SVLEN=', aster$SVLEN, ':NTLEN=', aster$NTLEN)
      aster = aster[ ,1:11]
      aster$cDNA_Change = ''
      aster$Refseq_mRNA_Id = ''
      aster$Codon_Change = ''
      aster$Protein_Change = ''
      apple = rbind(apple, aster)
      }
    }
  }

####################################################################################################################
# ScanITD:
aster.file = paste0(lock.dir, '/ScanITD/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  aster = aster[aster$N_Alt>thresh.n.alt, ]
  if (nrow(aster)>0)  {
    aster = cbind('ScanITD', aster)
    names(aster)[1] = 'Caller'
    aster$cDNA_Change = ''
    aster$Refseq_mRNA_Id = ''
    aster$Codon_Change = ''
    aster$Protein_Change = ''
    apple = rbind(apple, aster)
    }
  }

####################################################################################################################
# getITD:
aster.file = paste0(lock.dir, '/getITD/', sam, '.txt')
if (file.exists(aster.file))  {
  aster = read.table(aster.file, header=TRUE, quote='', sep='\t')
  aster = aster[aster$N_Alt>thresh.n.alt, ]
  if (nrow(aster)>0)  {
    aster = cbind('getITD', aster)
    names(aster)[1] = 'Caller'
    aster$cDNA_Change = ''
    aster$Refseq_mRNA_Id = ''
    aster$Codon_Change = ''
    aster$Protein_Change = ''
    apple = rbind(apple, aster)
    }
  }

####################################################################################################################
# Save the results and bat file to plot IGV for the sample if the ITD list is not empty.
if (length(apple)>0)  {
  apple = apple[order(apple$Batch, apple$Sample, apple$Start, -apple$N_Alt), ]
  apple$AF = round(apple$AF, 4)
  apple = apple[ ,c('Batch','Sample','Caller','Chrom','Start','End','Length','N_Alt','DP','AF','Type',
                    'cDNA_Change','Refseq_mRNA_Id','Codon_Change','Protein_Change')]

  # Save the results.
  write.table(apple, paste0(res.dir, '/Table/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

  # Create .bat file to plot zoom in IGV.
  apple = apple[ ,1:11]
  apple$start = round(apple$Start/35)
  apple$end = round(apple$End/35)
  apple = apple[ ,c(2,4:6,12:13)]
  idx = which(duplicated(apple[ ,c(1:2,5:6)]))
  if (length(idx)>0)  apple = apple[-idx, ]
  for (i in 1:nrow(apple))  {
    bam.dir = sub('ITD', 'BAM', lock.dir)
    bat = 'new'
    bat = rbind(bat, 'genome hg38')
    bat = rbind(bat, paste0('snapshotDirectory ', res.dir, '/IGV'))
    bat = rbind(bat, paste0('load ', bam.dir, '/', sam, '.bam'))
    bat = rbind(bat, 'maxPanelHeight 20000')
    bat = rbind(bat, 'preference SAM.SHOW_SOFT_CLIPPED TRUE')
    bat = rbind(bat, paste0('goto ', apple$Chrom[i], ':', apple$Start[i]-50, '-', apple$End[i]+50))
    bat = rbind(bat, paste0('snapshot ', sam, '_', apple$Start[i], '_', apple$End[i], '.png'))
    bat = rbind(bat, 'exit')
    write.table(bat, paste0(res.dir, '/IGV/', sam,'_', apple$Start[i], '_', apple$End[i], '.bat'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    }
  }  else  {
  apple = c('Batch', 'Sample', 'Caller', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', 'Type',
            'cDNA_Change', 'Refseq_mRNA_Id', 'Codon_Change', 'Protein_Change')
  write.table(t(apple), paste0(res.dir, '/Table/', sam, '.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
