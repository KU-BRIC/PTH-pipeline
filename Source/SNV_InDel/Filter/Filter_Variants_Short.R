####################################################################################################################
####################################################################################################################
# Filter variants - Short.
# Author: Haiying Kong
# Last Modified: 22 March 2021
####################################################################################################################
####################################################################################################################
options(stringsAsFactors=FALSE)
rm(list=ls())

library(xlsx)

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
res.dir = args[1]
sam = args[2]

####################################################################################################################
# Set values.
filter.style = 'Short'
thresh.n.alt = 4
thresh.maf = 0.01
var.classes = c('DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
                'Splice_Site', 'START_CODON_SNP', 'Translation_Start_Site')

####################################################################################################################
####################################################################################################################
# Read in the table with all variants.
var = read.table(paste0(res.dir, '/AllVariants/Callers_Wide/', sam, '.maf'), header=TRUE, quote='', sep='\t')

####################################################################################################################
####################################################################################################################
# Perform filtering.
####################################################################################################################
# Filter out technical errors.
####################################################################################################################
thresh.high.dp = read.table('/home/projects/cu_10184/projects/PTH/Reference/Filtering/Thresh_HighDP.txt',
                            header=FALSE, sep='\t')[1,1]
var = var[(var$t_alt_count>=thresh.n.alt & var$DP<=thresh.high.dp), ]

####################################################################################################################
# Filter out unwanted variant classes.
####################################################################################################################
var = var[((var$Chromosome %in% paste0('chr', c(1:22, 'X', 'Y'))) & var$Hugo_Symbol!='Unknown'), ]
var = var[var$Variant_Classification %in% var.classes, ]

####################################################################################################################
# Filter out polymorphysm with public database.
####################################################################################################################
# Identify variants for exclusion candidate.

# dbSNP:
idx.dbsnp.common = grep('1', var$dbSNP_COMMON)
idx.dbsnp.g5 = grep('true', var$dbSNP_G5)
idx.dbsnp.cfl = grep('true', var$dbSNP_CFL)
idx.dbsnp.gno = which(var$dbSNP_GNO=='true')
idx.dbsnp = unique(c(idx.dbsnp.common, idx.dbsnp.g5, idx.dbsnp.cfl, idx.dbsnp.gno))
# idx.dbsnp = unique(c(idx.dbsnp.common, idx.dbsnp.g5, idx.dbsnp.cfl))

# ExAC:
idx.exac = which(var$ClinVar_VCF_AF_EXAC > thresh.maf)
# ESP:
idx.esp = which(var$ClinVar_VCF_AF_ESP > thresh.maf)
# ClinVar_VCF_AF_TGP:
idx.tgp = which(var$ClinVar_VCF_AF_TGP > thresh.maf)
idx.maf = unique(c(idx.exac, idx.esp, idx.tgp))

# 1000G:
idx.1000g = grep('by1000genomes', var$dbSNP_Val_Status)

# ClinVar:
idx.clinvar.ex = unique(c(grep('Benign', var$ClinVar_VCF_CLNSIG), grep('Likely_benign', var$ClinVar_VCF_CLNSIG)))

# COSMIC:
idx.cosmic.ex = which(var$COSMIC_overlapping_mutations=='')

# idx.exclude = unique(c(idx.dbsnp.common, idx.dbsnp.g5, idx.dbsnp.cfl, idx.exac, idx.esp, idx.tgp, idx.1000g, idx.clinvar.ex, idx.cosmic.ex))
idx.exclude = unique(c(idx.dbsnp, idx.maf, idx.clinvar.ex, idx.cosmic.ex))

####################################################################################################################
# Identify variants for inclusion with highest priority.

# ClinVar:
idx.clinvar.in = unique(c(grep('Pathogenic', var$ClinVar_VCF_CLNSIG), grep('Likely_pathogenic', var$ClinVar_VCF_CLNSIG),
                        grep('risk_factor', var$ClinVar_VCF_CLNSIG), grep('drug_response', var$ClinVar_VCF_CLNSIG)))

# CIViC:               
civic = read.table('/home/projects/cu_10184/people/haikon/Reference/CIViC/hg38/01-Aug-2020-civic_accepted_and_submitted.vcf', header=FALSE, quote='', sep='\t')[ ,c(1,2,4,5)]
civic$Label = paste(civic$V1, as.character(civic$V2), civic$V4, civic$V5, sep='_')
idx.civic = which(var$Label %in% civic$Label)

idx.include = unique(c(idx.clinvar.in, idx.civic))

# Create one column to tag inclusion.
var$Inclusion = 0
var$Inclusion[idx.include] = 1

####################################################################################################################
# Combine inclusion and exclusion list with inclusion higher priority.

# Find index of variants for final exclusion.
idx = idx.exclude[!(idx.exclude %in% idx.include)]

# Update variant list.
var = var[-idx, ]

####################################################################################################################
# Filter out polymorphysm with our normal samples.
pon = read.table(paste0('/home/projects/cu_10184/projects/PTH/Reference/Filtering/', filter.style, '/Polymorphisms_PoN.txt'),
                 header=TRUE, sep='\t')
flag.var = paste(var$Hugo_Symbol, var$Chromosome, var$Start_Position, var$End_Position, var$Reference_Allele, var$Tumor_Seq_Allele2, sep='_')
flag.pon = apply(pon, 1, function(x) paste(x, collapse='_'))
idx = which(flag.var %in% flag.pon)

if (length(idx) > 0)  {
  var = var[-idx, ]
  flag.var = flag.var[-idx]
  }

####################################################################################################################
# Filter out region specific technical errors.
aster = read.table(paste0('/home/projects/cu_10184/projects/PTH/Reference/Filtering/', filter.style, '/RegionSpecificTechnicalError.txt'),
                   header=TRUE, sep='\t')
aster$flag = apply(aster[ ,1:6], 1, function(x) paste(x, collapse='_'))
idx.var = which(flag.var %in% aster$flag)

if (length(idx.var) > 0)  {
  i.del = c()
  idx.aster = match(flag.var[idx.var], aster$flag)
  for (i in 1:length(idx.var))  {
    if (abs(aster$UpperCenter[idx.aster[i]]-var$AF[idx.var[i]]) > abs(aster$LowerCenter[idx.aster[i]]-var$AF[idx.var[i]]))
      i.del = c(i.del, idx.var[i])
    }
  if (length(i.del) > 0)  var = var[-i.del, ]
  }

####################################################################################################################
# Save the list of variants after filtering.
write.table(var, paste0(res.dir, '/Filtered/', filter.style, '/', sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.xlsx(var, paste0(res.dir, '/Filtered/', filter.style, '/', sam, '.xlsx'), sheetName=paste0('Filtered_', filter.style), 
           row.names=FALSE, col.names=TRUE, append=FALSE)

####################################################################################################################
####################################################################################################################
