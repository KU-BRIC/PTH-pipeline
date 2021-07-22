####################################################################################################################
####################################################################################################################
# Get a list of normal samples sequenced with panel2.
# Author: Haiying Kong
# Last Modified: 21 July 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH/Meta')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Read in panel information, and get list of batches for panel2.
panel = read.table('BatchInfo.txt', header=TRUE, sep='\t')
batches = panel$Batch[panel$Target=='all_target_segments_covered_by_probes_Schmidt_Myeloid_TE-98545653_hg38_190919225222_updated-to-v2']

# Read in sample information, and get list of normal samples for panel2.
samples = read.table('SampleInfo.txt', header=TRUE, sep='\t')
normals = samples[((samples$Batch %in% batches) & (samples$Group=='NORMAL')), ]

# Save the list of variants after filtering.
write.table(normals, 'SampleInfo_Normal_Panel2.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
