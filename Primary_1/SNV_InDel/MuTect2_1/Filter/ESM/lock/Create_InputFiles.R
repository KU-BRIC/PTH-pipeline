####################################################################################################################
####################################################################################################################
# .
# Author: Haiying Kong
# Last Modified: 22 August 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get argument values from command line.
# args = commandArgs(trailingOnly=TRUE)
# batch_dir = args[1]
# sam = args[2]

dir.name = 'PTH'
batch = 'Primary_001'
sam = 'Horizon-CMP001'
temp.dir = paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/ESM/temp')


####################################################################################################################
####################################################################################################################
# Read in maf file.
maf = read.table(paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/TechError/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')
maf = maf[maf$Variant_Type=='SNP', ]
maf.tag = apply(maf[ ,c('Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')], 1, function(x) { x[2] = str_replace_all(as.character(x[2]), ' ', ''); paste(as.character(x), collapse='_')})
names(maf.tag) = NULL

# Read in vcf file.
vcf = read.table(paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Lock/SNV_InDel/MuTect2_1/vcf/', sam, '.vcf'), header=FALSE, quote='', sep='\t')
vcf.tag = apply(vcf[ ,c(1,2,4,5)], 1, function(x) { x[2] = str_replace_all(as.character(x[2]), ' ', ''); paste(x, collapse='_')})

# Create and save a new vcf file for the variats to process.
idx = match(maf.tag, vcf.tag)
vcf = vcf[idx, ]
write.table('##fileformat=VCFv4.2', paste0(temp.dir, '/vcf.vcf'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
names(vcf) = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SampleName')
suppressWarnings(write.table(vcf, paste0(temp.dir, '/vcf.vcf'), row.names=FALSE, col.names=TRUE, quote=FALSE, append=TRUE, sep='\t'))


vcf = readVcf(paste0(temp.dir, '/vcf.vcf'), "hg38")
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
coding = predictCoding(vcf, txdb, seqSource=Hsapiens)

maf$mutant = 

genes = unique(maf[ ,c('Hugo_Symbol', 'Chromosome')])

for (i in 1:nrow(genes))   {
  aster = maf[(maf$Hugo_Symbol==genes$Hugo_Symbol[i] & maf$Chromosome==genes$Chromosome[i]), ]
  protein.seq  = 

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)



####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

####################################################################################################################
write.table(long, paste0(res.dir, '/', sam, '_Long.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
if (nrow(short)>0)
  write.table(short, paste0(res.dir, '/', sam, '_Short.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
