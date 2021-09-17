####################################################################################################################
####################################################################################################################
# Run ESM for one sample.
# Author: Haiying Kong
# Last Modified: 28 August 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

options(stringsAsFactors=FALSE)
rm(list=ls())

library(VariantAnnotation)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
hasProteinData(edb)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
library(BSgenome.Hsapiens.UCSC.hg38)

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
dir.name = args[1]
batch = args[2]
sam = args[3]

esm.dir = paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript')

####################################################################################################################
####################################################################################################################
# Define function to check single amino acid substitution.
amino.acids = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
check.saas = function(p.change)  {
               prt = sub('p.', '', p.change)
               ref = substr(prt, 1, 1)
               alt = substr(prt, nchar(prt), nchar(prt))
               loc = substr(prt, 2, nchar(prt)-1)
               flag = as.integer((ref %in% amino.acids) & (alt %in% amino.acids) & !is.na(suppressWarnings(as.integer(loc))))
               return(list(flag, prt, ref, alt, loc))
               }

####################################################################################################################
####################################################################################################################
# Read in maf file.
maf = read.table(paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/TechError/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')

# Exclude variants with empty Protein_Change.
idx = which(maf$Protein_Change!='' & maf$Reference_Allele!='-' & maf$Tumor_Seq_Allele2!='-')

# Save non-SAAS.
write.table(maf[-idx, ], paste0(esm.dir, '/Original/', sam, '_rest.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# Candidate SAAS:
maf = maf[idx, ]

####################################################################################################################
# Create a vcf file for the amf file.
vcf.granges = maf[ ,c('Chromosome', 'Start_Position', 'End_Position', 'Strand')]
names(vcf.granges) = c('chr', 'start', 'end', 'strand')
vcf.granges = makeGRangesFromDataFrame(vcf.granges)

vcf.coldata = data.frame(Samples = 1)
row.names(vcf.coldata) = sam

vcf.header = scanVcfHeader(paste0('/home/projects/cu_10184/projects/', dir.name, '/BatchWork_1/', batch, '/Lock/SNV_InDel/MuTect2_1/vcf/', sam, '.vcf'))

vcf.fixed = DataFrame(REF = DNAStringSet(maf$Reference_Allele),
                      ALT = DNAStringSetList(maf$Tumor_Seq_Allele2))

vcf.info = maf[ ,193:217]

vcf = VCF(rowRanges=vcf.granges, colData=DataFrame(vcf.coldata), exptData=list(header=vcf.header), fixed=vcf.fixed, info=DataFrame(vcf.info))
####################################################################################################################

granges = maf[ ,c('Chromosome', 'Start_Position', 'End_Position', 'Strand')]
names(granges) = c('chr', 'start', 'end', 'strand')
granges = makeGRangesFromDataFrame(granges)

allele = DNAStringSet(maf$Tumor_Seq_Allele2)
coding = predictCoding(query=granges, subject=txdb, seqSource=Hsapiens, varAllele=allele)

coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)


  write.table(seq, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_seq.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  write.table(left, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_offset.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

  # maf file:
  write.table(maf[i, ], paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_acorn.maf'), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

  # input file.
  aster = as.data.frame(Protein_Change[i])
  names(aster) = 'Protein_Change'
  write.table(aster, paste0(esm.dir, '/temp/', sam, '_', as.character(i), '_input.csv'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=',')

  # Run ESM:
  system(paste0('sh /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2_1/Filter/ESM/OneTranscript/esm_one_variant.sh -d ', dir.name, ' -b ', batch, ' -s ', sam, ' -i ', as.character(i)))
  }

####################################################################################################################
####################################################################################################################
