####################################################################################################################
####################################################################################################################
# Combine all variants identified from VarDict, SNVer, and LoFreq.
# Author: Haiying Kong
# Last Modified: 29 March 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
lock.dir = args[1]
res.dir = args[2]
batch = args[3]
sam = args[4]

library(reshape2)

####################################################################################################################
####################################################################################################################
# Read in maf column names.
maf.cols = read.table('/home/projects/cu_10184/projects/PTH/Reference/MAF_Columns/MAF_Cols_trim0.txt', header=FALSE, sep='\t')[ ,1]

####################################################################################################################
####################################################################################################################
# Define function to modify variant classification for IDH2 missense mutation variants.
####################################################################################################################
IDH2_Update = function(maf)  {
  idx = which(maf$Genome_Change=='g.chr15:90088702C>T') 
  if (length(idx) > 0)  {
    maf$Variant_Classification[idx] = 'Missense_Mutation'
    maf$cDNA_Change[idx] = 'c.419G>A'
    maf$Protein_Change[idx] = 'p.R140Q'
    }
  idx = which(maf$Genome_Change=='g.chr15:90088703G>A')
  if (length(idx) > 0)  {
    maf$Variant_Classification[idx] = 'Missense_Mutation'
    maf$cDNA_Change[idx] = 'c.418C>T'
    maf$Protein_Change[idx] = 'p.R140W'
    }
  idx = which(maf$Genome_Change=='g.chr15:90088606C>T')
  if (length(idx) > 0)  {
    maf$Variant_Classification[idx] = 'Missense_Mutation'
    maf$cDNA_Change[idx] = 'c.515G>A'
    maf$Protein_Change[idx] = 'p.R172K'
    }
  return(maf)
  }

####################################################################################################################
####################################################################################################################
# Create large maf file by collecting all variants from all 3 callers.
####################################################################################################################
apple = c()
####################################################
# VarDict:
vcf = read.table(paste0(lock.dir, '/VarDict/vcf/', sam, '.vcf'), header=FALSE, sep='\t')[ ,c(2,4,5)]
maf = read.table(paste0(lock.dir, '/VarDict/maf/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols[-c(1:3,12:14)]]
maf$Batch = batch
maf$Sample = sam
maf$VarCall = 'VarDict'
maf$Start_vcf = vcf[ ,1]
maf$Ref_vcf = vcf[ ,2]
maf$Alt_vcf = vcf[ ,3]

# Update IDH2 variant information.
maf = IDH2_Update(maf)

# Update apple if maf is not empty.
if (nrow(maf)>0)  apple = rbind(apple, unique(maf[ ,maf.cols]))

####################################################
# SNVer:

# SNV:
vcf = read.table(paste0(lock.dir, '/SNVer/vcf/', sam, '.filter.vcf'), header=FALSE, sep='\t')[ ,c(2,4,5)]
maf = read.table(paste0(lock.dir, '/SNVer/maf/', sam, '.snv.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,c(maf.cols[-c(1:3,12:16)], 'AC', 'DP')]
maf$Batch = batch
maf$Sample = sam
maf$VarCall = 'SNVer'
maf$Start_vcf = vcf[ ,1]
maf$Ref_vcf = vcf[ ,2]
maf$Alt_vcf = vcf[ ,3]
maf$AF = maf$AC / maf$DP
maf$t_alt_count = maf$AF * maf$DP

# Update IDH2 variant information.
maf = IDH2_Update(maf)

# Update apple if maf is not empty.
if (nrow(maf)>0)  apple = rbind(apple, unique(maf[ ,maf.cols]))

# InDel:
vcf = read.table(paste0(lock.dir, '/SNVer/vcf/', sam, '.indel.filter.vcf'), header=FALSE, sep='\t')[ ,c(2,4,5)]
maf = read.table(paste0(lock.dir, '/SNVer/maf/', sam, '.indel.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,c(maf.cols[-c(1:3,12:16)], 'AC', 'DP')]
maf$Batch = batch
maf$Sample = sam
maf$VarCall = 'SNVer'
maf$Start_vcf = vcf[ ,1]
maf$Ref_vcf = vcf[ ,2]
maf$Alt_vcf = vcf[ ,3]
maf$AF = maf$AC / maf$DP
maf$t_alt_count = maf$AF * maf$DP

# Update IDH2 variant information.
maf = IDH2_Update(maf)

# Update apple if maf is not empty.
if (nrow(maf)>0)  apple = rbind(apple, unique(maf[ ,maf.cols]))

####################################################
# LoFreq:
vcf = read.table(paste0(lock.dir, '/LoFreq/vcf/', sam, '.vcf'), header=FALSE, sep='\t')[ ,c(2,4,5)]
maf = read.table(paste0(lock.dir, '/LoFreq/maf/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols[-c(1:3,12:14)]]
maf$Batch = batch
maf$Sample = sam
maf$VarCall = 'LoFreq'
maf$Start_vcf = vcf[ ,1]
maf$Ref_vcf = vcf[ ,2]
maf$Alt_vcf = vcf[ ,3]
maf$t_alt_count = maf$AF * maf$DP

# Update IDH2 variant information.
maf = IDH2_Update(maf)

# Update apple if maf is not empty.
if (nrow(maf)>0)  apple = rbind(apple, unique(maf[ ,maf.cols]))

####################################################################################################################
# Sort the result table.
apple = unique(apple)
apple = apple[order(apple$Batch, apple$Sample, apple$Hugo_Symbol,
                    apple$Chromosome, apple$Start_Position, apple$End_Position, apple$Strand,
                    apple$Reference_Allele, apple$Tumor_Seq_Allele1, apple$Tumor_Seq_Allele2, -apple$AF), ]

# Save the large table with all variants from all calling methods and all samples.
write.table(apple, paste0(res.dir, '/Callers_Long/', sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
# Reshape the table to wide form with one variant taking one row.
####################################################################################################################
# Identify rows with largest AF if multiple callings.
id.cols = names(apple)[-match(c('VarCall', 'AF', 'DP'), names(apple))]
id.cols.ess = c('Batch', 'Sample', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Strand',
                'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2')
idx = which(!duplicated(apple[ ,id.cols.ess]))

# Main rows in the wide table:
acorn = apple[idx, ]

####################################################################################################################
# Reshape for the first 12 columns in the long table to wide.
lace = apple[ ,c(id.cols.ess, 'VarCall', 'AF')]
# rm(list='apple')

formul = as.formula(paste0(paste(id.cols.ess, collapse='+'), ' ~ VarCall'))
lace = dcast(lace, formul, value.var='AF')
names(lace)[11:13] = paste('AF', names(lace)[11:13], sep='_')
lace[ ,11:13][is.na(lace[ ,11:13])] = 0

####################################################################################################################
# Merge lace and acorn.
acorn = merge(lace, acorn, by=id.cols.ess)
write.table(acorn, paste0(res.dir, '/Callers_Wide/', sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
