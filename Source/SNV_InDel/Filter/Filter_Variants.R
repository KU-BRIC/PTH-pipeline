####################################################################################################################
####################################################################################################################
# Filter variants.
# Author: Haiying Kong
# Last Modified: 19 December 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

options(stringsAsFactors=FALSE)
rm(list=ls())

library(parallel)
library(xlsx)

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
res.dir = args[1]
sam = args[2]
scheme_dir = args[3]
scheme_name = args[4]

####################################################################################################################
# Set values.
var.classes = c('DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
                'Splice_Site', 'START_CODON_SNP', 'Translation_Start_Site')

####################################################################################################################
####################################################################################################################
# Read in data for the filtering scheme.
ref.dir = paste0('/home/projects/cu_10184/projects/', scheme_dir, '/Reference/Filtering/', scheme_name)
thresh = read.table(paste0(ref.dir, '/Thresholds.txt'), header=TRUE, sep='\t')
pon = read.table(paste0(ref.dir, '/SNPs_PoN.txt'), header=TRUE, sep='\t')
error = read.table(paste0(ref.dir, '/RegionSpecificTechnicalError.txt'), header=TRUE, sep='\t')

####################################################################################################################
####################################################################################################################
# Read in the table with all variants.
var = read.table(paste0(res.dir, '/AllVariants/Callers_Wide/', sam, '.maf'), header=TRUE, quote='', sep='\t')

####################################################################################################################
####################################################################################################################
# thresh_dp_high_med_fold = 10
# thresh_dp_high = quantile(var$DP,probs=0.5) * thresh_dp_high_med_fold

thresh_dp_high_IQR_fold = 2
thresh_dp_high = quantile(var$DP,probs=0.75) + IQR(var$DP) * thresh_dp_high_IQR_fold

####################################################################################################################
####################################################################################################################
# Perform filtering.
####################################################################################################################
# Filter out technical errors.
####################################################################################################################
var = var[(var$t_alt_count>=thresh$thresh_n_alt & var$DP>=thresh$thresh_dp_low & var$DP<=thresh_dp_high), ]

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
idx.exac = which(var$ClinVar_VCF_AF_EXAC > thresh$thresh_maf_db)
# ESP:
idx.esp = which(var$ClinVar_VCF_AF_ESP > thresh$thresh_maf_db)
# ClinVar_VCF_AF_TGP:
idx.tgp = which(var$ClinVar_VCF_AF_TGP > thresh$thresh_maf_db)
idx.maf = unique(c(idx.exac, idx.esp, idx.tgp))

# 1000G:
idx.1000g = grep('by1000genomes', var$dbSNP_Val_Status)

# ClinVar:
idx.clinvar.ex = unique(c(grep('Benign', var$ClinVar_VCF_CLNSIG), grep('Likely_benign', var$ClinVar_VCF_CLNSIG)))

# COSMIC:
if (scheme_name=='Short')  {
  idx.cosmic.ex = which(var$COSMIC_overlapping_mutations=='')
  }  else  {
  idx.cosmic.ex = which(var$COSMIC_tissue_types_affected=='')
  }

###################################
# Polymorphisms identified from our normal samples.
var.flag = paste(var$Chromosome, as.character(var$Start_Position), as.character(var$End_Position), var$Reference_Allele, var$Tumor_Seq_Allele2, sep='_')
flag.pon = paste(pon$Chrom, as.character(pon$Start), as.character(pon$End), pon$Ref, pon$Alt, sep='_')
idx.pon = which(var.flag %in% flag.pon)

###################################
# Black list:
bed = read.table('/home/projects/cu_10184/projects/PTH/Reference/Filtering/BlackList/black.bed', header=TRUE, sep='\t')
bed = bed[(bed$Sample=='Global' | bed$Sample==sam), ]
n.cores = detectCores()
ans1 = mclapply(1:nrow(var), function(i)  {
                               idx = which(bed$Chrom==var$Chromosome[i] & bed$Start<=var$Start_Position[i] & bed$End>=var$End_Position[i])
                               if (length(idx)>0)
                                 flag = 1   else
                                 flag = 0   
                               flag
                               },
                mc.cores = n.cores)
ans1 = which(ans1==1)

vcf = read.table('/home/projects/cu_10184/projects/PTH/Reference/Filtering/BlackList/black.vcf', header=TRUE, sep='\t')
vcf = vcf[(vcf$Sample=='Global' | vcf$Sample==sam), ]
flag.vcf = apply(vcf, 1, function(x) paste(x, collapse='_'))
var.flag1 = paste(var$Chromosome, var$Start_Position, var$Reference_Allele, var$Tumor_Seq_Allele2, sep='_')
ans2 = which(var.flag1 %in% flag.vcf)

idx.black = c(ans1, ans2)

###################################
# Find index list for exclusion.
idx.exclude = unique(c(idx.dbsnp, idx.maf, idx.clinvar.ex, idx.cosmic.ex, idx.pon, idx.black))

####################################################################################################################
# Identify variants for inclusion with highest priority.

# ClinVar:
idx.clinvar.in = unique(c(grep('Pathogenic', var$ClinVar_VCF_CLNSIG), grep('Likely_pathogenic', var$ClinVar_VCF_CLNSIG),
                        grep('risk_factor', var$ClinVar_VCF_CLNSIG)))

# CIViC:               
civic = read.table('/home/projects/cu_10184/people/haikon/Reference/CIViC/hg38/01-Aug-2020-civic_accepted_and_submitted.vcf', header=FALSE, quote='', sep='\t')[ ,c(1,2,4,5)]
civic$flag = paste(civic$V1, as.character(civic$V2), civic$V4, civic$V5, sep='_')
idx.civic = which(var.flag %in% civic$flag)

###################################
# White list:
bed = read.table('/home/projects/cu_10184/projects/PTH/Reference/Filtering/WhiteList/white.bed', header=TRUE, sep='\t')
bed = bed[(bed$Sample=='Global' | bed$Sample==sam), ]
n.cores = detectCores()
ans1 = mclapply(1:nrow(var), function(i)  {
                               idx = which(bed$Chrom==var$Chromosome[i] & bed$Start<=var$Start_Position[i] & bed$End>=var$End_Position[i])
                               if (length(idx)>0)
                                 flag = 1   else
                                 flag = 0
                               flag
                               },
                mc.cores = n.cores)
ans1 = which(ans1==1)

vcf = read.table('/home/projects/cu_10184/projects/PTH/Reference/Filtering/WhiteList/white.vcf', header=TRUE, sep='\t')
vcf = vcf[(vcf$Sample=='Global' | vcf$Sample==sam), ]
flag.vcf = apply(vcf, 1, function(x) paste(x, collapse='_'))
var.flag1 = paste(var$Chromosome, var$Start_Position, var$Reference_Allele, var$Tumor_Seq_Allele2, sep='_')
ans2 = which(var.flag1 %in% flag.vcf)

idx.white = c(ans1, ans2)

idx.include = unique(c(idx.clinvar.in, idx.civic, idx.white))

####################################################################################################################
# Combine inclusion and exclusion list with inclusion higher priority.

# Find index of variants for final exclusion.
idx = idx.exclude[!(idx.exclude %in% idx.include)]

# Update variant list.
var = var[-idx, ]

####################################################################################################################
# Filter out region specific technical errors.
error$flag = apply(error[ ,1:6], 1, function(x) paste(x, collapse='_'))
idx.var = which(var.flag %in% error$flag)

if (length(idx.var) > 0)  {
  i.del = c()
  idx.error = match(var.flag[idx.var], error$flag)
  for (i in 1:length(idx.var))  {
    if (abs(error$UpperCenter[idx.error[i]]-var$AF[idx.var[i]]) > abs(error$LowerCenter[idx.error[i]]-var$AF[idx.var[i]]))
      i.del = c(i.del, idx.var[i])
    }
  if (length(i.del) > 0)  var = var[-i.del, ]
  }

var = var[ ,-ncol(var)]

####################################################################################################################
# Save the list of variants after filtering.
write.table(var, paste0(res.dir, '/Filtered/', scheme_name, '/', sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
if (nrow(var)>0)  {
  write.xlsx(var, paste0(res.dir, '/Filtered/', scheme_name, '/', sam, '.xlsx'), sheetName=paste0('Filtered_', scheme_name),
             row.names=FALSE, col.names=TRUE, append=FALSE)
  }

####################################################################################################################
####################################################################################################################
