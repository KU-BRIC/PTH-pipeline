####################################################################################################################
####################################################################################################################
# Create amplicon_kayser.tsv an annotation file for getITD.
# Author: Haiying Kong
# Last Modified: 15 April 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

####################################################################################################################
####################################################################################################################
# Load exon database, create target file and filter for exons on target.
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
exon = exons(txdb)
target = GRanges(seqnames = rep('chr13', 3),
                 ranges = IRanges(start=c(28033736,28033931,28034150), end=c(28034141,28034364,28034557)))
exon = as.data.frame(subsetByOverlaps(exon, target))
target = as.data.frame(target)

####################################################################################################################
# Create a data.frame for all positions on target as rows.
aster = (target$end[nrow(target)]:target$start[1])
apple = data.frame(amplicon_bp = 1:length(aster),
                   region = NA,
                   chr13_bp = aster,
                   transcript_bp = '.',
                   protein_as = '.')

####################################################################################################################
# Annotate for 3 exons.
i.exon = 15
for (i in nrow(exon):1)  {
  ii = which(apple$chr13_bp>=exon$start[i] & apple$chr13_bp<=exon$end[i])
  apple$region[ii] = paste0('exon', i.exon)
  i.exon = i.exon - 1
  apple$transcript_bp[ii] = as.character(1:length(ii))
  apple$protein_as[ii] = as.character(as.integer((as.integer(apple$transcript_bp[ii])-0.1)/3) + 1)
  }
  
apple$region[apple$chr13_bp>exon$end[3]] = 'intron_15'
apple$region[apple$chr13_bp>exon$end[2] & apple$chr13_bp<exon$start[3]] = 'intron_14'
apple$region[apple$chr13_bp>exon$end[1] & apple$chr13_bp<exon$start[2]] = 'intron_13'
apple$region[apple$chr13_bp<exon$start[1]] = 'intron_12'

# Add one more row after one last row.
n = nrow(apple)
apple[(n+1), ] = apple[n, ]
apple$chr13_bp[n+1] = apple$chr13_bp[n+1] - 1
apple$amplicon_bp[n+1] = apple$amplicon_bp[n+1] + 1

write.table(apple, 'Reference/ITD/FLT3_1/getITD/amplicon_kayser.tsv', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
