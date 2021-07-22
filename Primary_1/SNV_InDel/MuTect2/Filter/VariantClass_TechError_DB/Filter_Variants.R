####################################################################################################################
####################################################################################################################
# Filter variants.
# Author: Haiying Kong
# Last Modified: 21 July 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
batch_dir = args[1]
sam = args[2]

####################################################################################################################
# Set values.
####################################################################################################################
# Variant classes:
var.classes = c('DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
                'Splice_Site', 'START_CODON_SNP', 'Translation_Start_Site')

# Thresholds for technical errors:
thresh_n_alt = 4
thresh_dp_low = 200
thresh_dp_high = 2873

# Thresholds for public database:
thresh_maf_db = 0.01

####################################################################################################################
####################################################################################################################
# Read in the table with all variants.
var = read.table(paste0(batch_dir, '/Lock/SNV_InDel/MuTect2/maf/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')

if (nrow(var)==0)  quit(save='no', status=0)

####################################################################################################################
####################################################################################################################
# Perform filtering.
####################################################################################################################
# Filter variant classes of our interest.
####################################################################################################################
scheme_name = 'VariantClass'

var = var[((var$Chromosome %in% paste0('chr', c(1:22, 'X', 'Y'))) & var$Hugo_Symbol!='Unknown'), ]
var = var[var$Variant_Classification %in% var.classes, ]

if (nrow(var)==0)  quit(save='no', status=0)

res.dir = paste0(batch_dir, '/Result/SNV_InDel/MuTect2/', scheme_name)
write.table(var, paste0(res.dir, '/', sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
# Filter out technical errors.
####################################################################################################################
scheme_name = 'TechError'

var = var[(var$t_alt_count>=thresh_n_alt & var$DP>=thresh_dp_low & var$DP<=thresh_dp_high), ]

if (nrow(var)==0)  quit(save='no', status=0)

res.dir = paste0(batch_dir, '/Result/SNV_InDel/MuTect2/', scheme_name)
write.table(var, paste0(res.dir, '/', sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
# Filter out polymorphisms or non-pathogenic variants from public database.
####################################################################################################################
scheme_name = 'DB'

###################################
# Identify variants for exclusion candidate.

# dbSNP:
idx.dbsnp.common = grep('1', var$dbSNP_COMMON)
idx.dbsnp.g5 = grep('true', var$dbSNP_G5)
idx.dbsnp.cfl = grep('true', var$dbSNP_CFL)
idx.dbsnp.gno = which(var$dbSNP_GNO=='true')
idx.dbsnp = unique(c(idx.dbsnp.common, idx.dbsnp.g5, idx.dbsnp.cfl, idx.dbsnp.gno))
# idx.dbsnp = unique(c(idx.dbsnp.common, idx.dbsnp.g5, idx.dbsnp.cfl))

# ExAC:
idx.exac = which(var$ClinVar_VCF_AF_EXAC > thresh_maf_db)
# ESP:
idx.esp = which(var$ClinVar_VCF_AF_ESP > thresh_maf_db)
# ClinVar_VCF_AF_TGP:
idx.tgp = which(var$ClinVar_VCF_AF_TGP > thresh_maf_db)
idx.maf = unique(c(idx.exac, idx.esp, idx.tgp))

# 1000G:
idx.1000g = grep('by1000genomes', var$dbSNP_Val_Status)

# ClinVar:
idx.clinvar.ex = unique(c(grep('Benign', var$ClinVar_VCF_CLNSIG), grep('Likely_benign', var$ClinVar_VCF_CLNSIG)))

# COSMIC:
idx.cosmic.ex.long = which(var$COSMIC_tissue_types_affected=='')
idx.cosmic.ex.short = which(var$COSMIC_overlapping_mutations=='')

# Find index list for exclusion.
idx.exclude.long = unique(c(idx.dbsnp, idx.maf, idx.clinvar.ex, idx.cosmic.ex.long))
idx.exclude.short = unique(c(idx.dbsnp, idx.maf, idx.clinvar.ex, idx.cosmic.ex.short))

###################################
# Identify variants for inclusion.

# ClinVar:
idx.clinvar.in = unique(c(grep('Pathogenic', var$ClinVar_VCF_CLNSIG), grep('Likely_pathogenic', var$ClinVar_VCF_CLNSIG),
                        grep('risk_factor', var$ClinVar_VCF_CLNSIG), grep('drug_response', var$ClinVar_VCF_CLNSIG)))

# CIViC:               
civic = read.table('/home/projects/cu_10184/people/haikon/Reference/CIViC/hg38/01-Aug-2020-civic_accepted_and_submitted.vcf', header=FALSE, quote='', sep='\t')[ ,c(1,2,4,5)]
civic$Label = paste(civic$V1, as.character(civic$V2), civic$V4, civic$V5, sep='_')
idx.civic = which(var$Label %in% civic$Label)

# Find index list for inclusion.
idx.include = unique(c(idx.clinvar.in, idx.civic))

###################################
# Combine inclusion and exclusion list with inclusion higher priority.

# Find index of variants for final exclusion.
idx.long = idx.exclude.long[!(idx.exclude.long %in% idx.include)]
idx.short = idx.exclude.short[!(idx.exclude.short %in% idx.include)]

###################################
# Filter variants.
if (length(idx.long)>0)
  long = var[-idx.long, ]  else
  long = var
if (length(idx.short)>0)
  short = var[-idx.short, ]  else
  short = var

res.dir = paste0(batch_dir, '/Result/SNV_InDel/MuTect2/', scheme_name)
if (nrow(long)>0)
  write.table(long, paste0(res.dir, '/', sam, '_Long.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
if (nrow(short)>0)
  write.table(short, paste0(res.dir, '/', sam, '_Short.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
