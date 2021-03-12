####################################################################################################################
####################################################################################################################
# Get coverage frequency and plot density for all samples.
# Author: Haiying Kong modified by Balthasar
# Last Modified: 4 March 2017
####################################################################################################################
####################################################################################################################
args = commandArgs(trailingOnly=TRUE)

bamdir=args[1]
sample=args[2]
bedfile=args[3]
output_dir=args[4]

setwd(output_dir)

library(HelloRanges)

####################################################################################################################
####################################################################################################################
# Run the pipeline for QC of sequence data.
bam.file = paste(bamdir,sample,sep="/")
genome = Seqinfo(genome = NA_character_)
gr_a = import(bedfile, genome=genome)
gr_b = import(bam.file, genome=genome)

cov = unname(coverage(gr_b)[gr_a])
all_cov = unlist(cov)
df = as.data.frame(table(Coverage=all_cov))
df$Fraction = df$Freq/length(all_cov)

write.table(df, paste("FreqTable/", sub(".bam", "", sample), ".txt", sep=""),
        row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

png(paste("DensityPlot/", sub(".bam", "", sample), ".png", sep=""))
plot(Fraction~as.integer(Coverage), df, type="s",
	xlab="depth of coverage", ylab = "probability",
	main=sub(".bam", "", sample))
dev.off()

####################################################################################################################
####################################################################################################################

