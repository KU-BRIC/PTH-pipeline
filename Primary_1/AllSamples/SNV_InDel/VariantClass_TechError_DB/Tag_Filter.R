####################################################################################################################
####################################################################################################################
# Tag apple with known 
# Author: Haiying Kong
# Last Modified: 26 August 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Read in known Horizon and REDCap variants.
tag = read.table('/home/projects/cu_10184/projects/PTH/AllBatches/Result/SNV_InDel/ClinicTag/CalledReferenceVariants.txt',
                 header=TRUE, sep='\t')
tag = tag[ ,c('Batch', 'Sample', 'Hugo_Symbol', 'Chrom', 'Start', 'Ref', 'Alt', 'Refseq_ID', 'cDNA_Change', 'Protein_Change', 'AF_Clinic', 'Note')]
pat = sapply(tag$Sample, function(x)  unlist(strsplit(x, '-'))[1])
pat = as.vector(pat)
tag.label = paste(tag$pat, tag$Protein_Change, sep='-')

Protein_Change = gsub('p.', '', tag$Protein_Change)
amino.acids = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
flag = sapply(Protein_Change, function(x)  {
                                ref = substr(x, 1, 1)
                                alt = substr(x, nchar(x), nchar(x))
                                loc = substr(x, 2, nchar(x)-1)
                                flag = as.integer((ref %in% amino.acids) & (alt %in% amino.acids) & !is.na(suppressWarnings(as.integer(loc))))
                                })
idx = which(as.vector(flag==1))
tag = tag[idx, ]
tag.label = paste(pat[idx], tag$Protein_Change, sep='-')

####################################################################################################################
####################################################################################################################
#
####################################################################################################################
# Read in table with all SAAS.
maf = read.table('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/Apple_SAAS.txt',
                 header=TRUE, quote='', fill=TRUE, sep='\t')
pat = sapply(maf$Sample, function(x)  unlist(strsplit(x, '-'))[1])
pat = as.vector(pat)
maf.label = paste(pat, maf$Protein_Change, sep='-')


rest = read.table('/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/Apple_rest.txt',
                  header=TRUE, quote='', fill=TRUE, sep='\t')


# Save the results.
write.table(apple, 'AllBatches/Result/SNV_InDel/Filtered/Variants_VariantClass_TechErrpr_PoN.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
