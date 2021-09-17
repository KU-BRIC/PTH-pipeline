####################################################################################################################
####################################################################################################################
# Run SomaticSignatures.
# Author: Haiying Kong
# Last Modified: 7 September 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(ggdendro)
library(sva)

####################################################################################################################
# Set parameters.
proj.name = 'PTH'

####################################################################################################################
####################################################################################################################
# Read in maf file of all samples.
maf.file = paste0('/home/projects/cu_10184/projects/', proj.name, '/AllBatches_1/Lock/SNV_InDel/MuTect2_1/Apple.maf')
maf = SomaticCancerAlterations:::.read_maf(maf.file)
nucleo = c('A', 'C', 'G', 'T')
maf = maf[((maf$Reference_Allele %in% nucleo) & (maf$Tumor_Seq_Allele2 %in% nucleo)), ]

# Create input object.
vr = VRanges(seqnames = Rle(maf$Chromosome),
             ranges = IRanges(start=maf$Start_Position, end=maf$End_Position),
             ref = maf$Reference_Allele,
             alt = maf$Tumor_Seq_Allele2,
             sampleNames = maf$Tumor_Sample_Barcode,
             seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38),
             study = rep('PTH', nrow(maf)))

# Run SomaticSignatures.
pdf(paste0('/home/projects/cu_10184/projects/', proj.name, '/AllBatches_1/Result/SNV_InDel/MuTect2_1/SummaryVisual/SomaticSignatures.pdf'))

motifs = mutationContext(vr, BSgenome.Hsapiens.UCSC.hg38)
head(motifs)
mm = motifMatrix(motifs, group="study", normalize=FALSE)

head(round(mm, 4))
plotMutationSpectrum(motifs, "study")

n_sigs = 5
sigs = identifySignatures(mm, n_sigs, nmfDecomposition)
sigs

n_sigs = 2:8

gof = assessNumberSignatures(mm, n_sigs, nReplicates = 5)
plotNumberSignatures(gof)
plotSignatureMap(sigs) + ggtitle("Somatic Signatures: NMF - Heatmap")
plotSignatures(sigs) + ggtitle("Somatic Signatures: NMF - Barchart")
plotObservedSpectrum(sigs)
plotFittedSpectrum(sigs)
plotSampleMap(sigs)
plotSamples(sigs)

p = plotSamples(sigs)
p = p + theme(legend.position = "none")
p = p + xlab("Studies")
p = p + ggtitle("Somatic Signatures in TGCA WES Data")
p = p + scale_fill_brewer(palette = "Blues")
p = p + theme(axis.text.x = element_text(size = 9))
clu_motif = clusterSpectrum(mm, "motif")
p = ggdendrogram(clu_motif, rotate = TRUE)


# sca_anno = as.data.frame(lapply(sca_metadata, unlist))
model_null = model.matrix(~ 1, sca_anno)
sca_mm_batch = ComBat(sca_mm, batch = sca_anno$Sequence_Source, mod = model_null)
k = 3
n = 1e4
hs_chrs = as(seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5), "GRanges")
hs_chrs = keepSeqlevels(hs_chrs, c(1:22, "X", "Y"), pruning.mode = "coarse")

k3_hs_chrs = kmerFrequency(BSgenome.Hsapiens.1000genomes.hs37d5, n, k, hs_chrs)
k3_hs_chrs

k = 3
n = 1e4

hs_exons = reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
hs_exons = ncbi(keepStandardChromosomes(hs_exons))

k3_exons = kmerFrequency(BSgenome.Hsapiens.1000genomes.hs37d5, n, k, hs_exons)
data(kmers)
norms = k3wg / k3we
head(norms)

dev.off()

####################################################################################################################
####################################################################################################################
