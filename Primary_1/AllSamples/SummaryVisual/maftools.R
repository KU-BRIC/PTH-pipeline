####################################################################################################################
####################################################################################################################
# Run maftools.
# Author: Haiying Kong
# Last Modified: 7 September 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(maftools)

####################################################################################################################
# Set parameters.
proj.name = 'PTH'
n.var = 50
####################################################################################################################
####################################################################################################################
# Clean old folder and create new one to save results.
res.dir = paste0('/home/projects/cu_10184/projects/', proj.name, '/AllBatches_1/Result/SNV_InDel/MuTect2_1/SummaryVisual/maftools')
if (dir.exists(res.dir))  unlink(res.dir, recursive=TRUE)
dir.create(res.dir)
dir.create(paste0(res.dir, '/SummaryStats'))
dir.create(paste0(res.dir, '/Plots'))

####################################################################################################################
# Read in maf file of all samples.
maf.file = paste0('/home/projects/cu_10184/projects/', proj.name, '/AllBatches_1/Lock/SNV_InDel/MuTect2_1/Apple_TechErrorFiltered.maf')
maf = read.maf(maf=maf.file)

####################################################################################################################
# Run maftools for summary stats.
getSampleSummary(maf)
getGeneSummary(maf)
getFields(maf)
write.mafSummary(maf=maf, basename=paste0(res.dir, '/SummaryStats/TechErrorFiltered'))

####################################################################################################################
# Run maftools for visualization.

# maf summary plots.
pdf(paste0(res.dir, '/Plots/Summary.pdf'))
plotmafSummary(maf=maf, rmOutlier=TRUE, addStat='median', dashboard=TRUE, titvRaw=FALSE)
dev.off()

# oncoplots for genes of our interest.
pdf(paste0(res.dir, '/Plots/oncoplots.pdf'))
top.var = read.table('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/TopVariants_SAAS_ESM_2_byVariant.maf',
                     header=TRUE, quote='', fill=TRUE, sep='\t')
top.genes = unique(top.var$Hugo_Symbol[1:n.var])
aster = subsetMaf(maf=maf, genes=top.genes)
oncoplot(maf=aster)
dev.off()

# transition and transversions.
pdf(paste0(res.dir, '/Plots/Transition_Transversion.pdf'))
maf.titv = titv(maf=maf, plot=FALSE, useSyn=TRUE)
plotTiTv(res=maf.titv)
dev.off()

# Lollipop plots for amino acid changes.
pdf(paste0(res.dir, '/Plots/Lollipop.pdf'))
for (top.gene in top.genes)  {
  label.pos = unique(aster@data$Protein_Change[(aster@data$Hugo_Symbo==top.gene & aster@data$Variant_Classification=='Missense_Mutation')])
  label.pos = suppressWarnings(as.numeric(substr(label.pos, 4, nchar(label.pos)-1)))
  label.pos = label.pos[!is.na(label.pos)]
  lollipopPlot(maf=aster, gene=top.gene, AACol='Protein_Change', showMutationRate=TRUE, labelPos=label.pos)
  }
dev.off()

# plotProtein(gene='TP53', refSeqID='NM_000546')

# Rainfall plots.
pdf(paste0(res.dir, '/Plots/Rainfall.pdf'))
samples = unique(maf@data$Tumor_Sample_Barcode)
for (sam in samples)  {
  rainfallPlot(maf=subsetMaf(maf=maf,tsb=sam), detectChangePoints=TRUE, pointSize=0.4)
  }
dev.off()

# Compare mutation load against TCGA cohorts.
pdf(paste0(res.dir, '/Plots/Compare_MutationLoad.pdf'))
maf.mutload = tcgaCompare(maf=maf, cohortName='PTH', logscale=TRUE, capture_size=50)
maf@data$AF = 100*maf@data$AF
plotVaf(maf=maf, vafCol='AF')
dev.off()

# Plot VAF.
pdf(paste0(res.dir, '/Plots/Genes_VAF.pdf'))
plotVaf(maf=subsetMaf(maf=maf, genes=top.genes), vafCol='AF')
dev.off()

####################################################################################################################
####################################################################################################################
