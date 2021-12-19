####################################################################################################################
####################################################################################################################
# Get summary metrics table from QC.
# Author: Haiying Kong
# Last Modified: 19 December 2021
####################################################################################################################
####################################################################################################################
options(stringsAsFactors=FALSE)
rm(list=ls())

# Wanted columns:
apple.cols = c()

# Get passed argument values.
args = commandArgs(trailingOnly=TRUE)
dir_name = args[1]
batch = args[2]
sam = args[3]

setwd(paste0('/home/projects/cu_10184/projects/', dir_name))

# Set parameters.
doc.dir = paste0('/home/projects/cu_10184/projects/', dir_name, '/BatchWork/', batch, '/Lock/DepthOfCoverage/FreqTable')
qc.dir = paste0('/home/projects/cu_10184/projects/', dir_name, '/QC/Result/FASTQuick/ByBatch/', batch, '/', sam)

####################################################################################################################
####################################################################################################################
# Read in information from DOC and FASTQuick output.
dp = read.table(paste0(doc.dir, '/', sam, '.txt'), header=TRUE, sep='\t')
seq.summ = read.table(paste0(qc.dir, '/', sam, '.Sequence.csv'), header=TRUE, sep=',')
summ = read.table(paste0(qc.dir, '/', sam, '.Summary'), header=TRUE, sep=':')

# Extract critical information.
apple = data.frame(ReadLength = as.integer(seq.summ$ReadLength[1]),
                   DP_Median = round(sum(dp$Coverage * dp$Freq) / sum(dp$Freq)),
                   DP_GT_100 = round(sum(dp$Freq[dp$Coverage > 100]) / sum(dp$Freq), 4),
                   DP_GT_500 = round(sum(dp$Freq[dp$Coverage > 500]) / sum(dp$Freq), 4),
                   DP_GT_900 = round(sum(dp$Freq[dp$Coverage > 900]) / sum(dp$Freq), 4),
                   Median_IS_DP_GT_500 = as.integer(summ[15,2]),
                   Median_IS_DP_GT_300 = as.integer(summ[16,2]),
                   Q20_Frac = round(as.numeric(summ[11,2]), 4),
                   Q30_Frac = round(as.numeric(summ[12,2]), 4),
                   PCR_Dup_Rate = round(as.numeric(unlist(strsplit(summ[2,2],'\\['))[1]), 4),
                   Contamination = round(as.numeric(summ[17,2]), 4)
                   )

####################################################################################################################
# Save the results.
write.table(apple, paste0(qc.dir, '/OurSummary.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
