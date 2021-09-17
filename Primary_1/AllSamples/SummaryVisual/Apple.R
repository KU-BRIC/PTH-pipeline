####################################################################################################################
####################################################################################################################
# Create trimmed maf file for all variants identified from all samples.
# Author: Haiying Kong
# Last Modified: 7 September 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

####################################################################################################################
# Thresholds for technical errors:
thresh_n_alt = 4
thresh_dp_low = 200
thresh_dp_high = 2873

# Set parameters.
proj.name = 'PTH'
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

####################################################################################################################
# Read in the names of columns that will be kept.
maf.cols = c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
             'Variant_Classification', 'Variant_Type', 'Protein_Change', 'AF', 'DP')

####################################################################################################################
####################################################################################################################
# Create large maf file by concatenating all trimmed maf files from all samples.
####################################################################################################################
# Collect all variants from all samples with technical error, variant class and PoN filtered.
apple = c()
apple.filtered = c()

for (batch in batches)  {
  maf.dir = paste0('/home/projects/cu_10184/projects/', proj.name, '/BatchWork_1/', batch, '/Lock/SNV_InDel/MuTect2_1/maf')
  maf.files = dir(maf.dir, pattern='.maf')
  for (maf.file in maf.files)   {
    # All variants:
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
    maf$Tumor_Sample_Barcode = sub('.maf', '', maf.file)
    apple = rbind(apple, maf[ ,maf.cols])

    # Technical errors filtered:
    maf = maf[(maf$t_alt_count>=thresh_n_alt & maf$DP>=thresh_dp_low & maf$DP<=thresh_dp_high), ]
    if (nrow(maf)>0)  apple.filtered = rbind(apple.filtered, maf[ ,maf.cols])
    }
  }

# Save the maf files.
file.name = paste0('/home/projects/cu_10184/projects/', proj.name, '/AllBatches_1/Lock/SNV_InDel/MuTect2_1/Apple.maf')
write.table(apple, file.name, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

file.name = paste0('/home/projects/cu_10184/projects/', proj.name, '/AllBatches_1/Lock/SNV_InDel/MuTect2_1/Apple_TechErrorFiltered.maf')
write.table(apple, file.name, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
