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

library(Biostrings)

####################################################################################################################
####################################################################################################################
# Read in the fasta file.
fa = read.table('Reference/ITD/FLT3_1/getITD/amplicon.txt', header=FALSE, sep='\t')
dna = fa[-grep('^>', fa[ ,1]), 1]
dna = paste(dna, collapse='')

# Reverse the string and convert to lower case because FLT3 is encoded on negative strand.
dna = strsplit(dna, NULL)[[1]]
dna = paste(dna, collapse='')
dna = DNAString(dna)
dna = reverseComplement(dna)
dna = as.character(dna)

# Save the results.
write.table(dna, 'Reference/ITD/FLT3_1/getITD/amplicon.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
