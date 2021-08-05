####################################################################################################################
####################################################################################################################
# Create large table for all variants identified from all samples, and filtered with variant class and technical error.
# Author: Haiying Kong
# Last Modified: 24 July 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

####################################################################################################################
# Set parameters.
proj.name = 'PTH'
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

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
# Read in the names of columns that will be kept.
# maf.cols = read.table('/home/projects/cu_10184/projects/PTH/Reference/MAF_Columns/MAF_Cols_1.txt', header=FALSE, sep='\t')[ ,1]

# Read in PoN.
pon.file = '/home/projects/cu_10184/projects/PTH/Reference/Filtering/Medium/SNPs_PoN.txt'
pon = read.table(pon.file, header=TRUE, sep='\t')
flag.pon = apply(pon, 1, function(x) paste(x, collapse='_'))

####################################################################################################################
####################################################################################################################
# Create large maf file by concatenating all maf files from all samples.
####################################################################################################################
# Collect all variants from all samples with technical error, variant class and PoN filtered.
apple = c()

for (batch in batches)  {
  maf.dir = paste0('BatchWork/', batch, '/Result/SNV_InDel/AllVariants/Callers_Wide')
  maf.files = dir(maf.dir, pattern='.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
    if (nrow(maf)>0)  {

      # Filter with variant classes.
      maf = maf[((maf$Chromosome %in% paste0('chr', c(1:22, 'X', 'Y'))) & maf$Hugo_Symbol!='Unknown'), ]
      maf = maf[maf$Variant_Classification %in% var.classes, ]

      # Filter technical errors.
      maf = maf[(maf$t_alt_count>=thresh_n_alt & maf$DP>=thresh_dp_low & maf$DP<=thresh_dp_high), ]

      # Filter with PoN.
      flag.maf = paste(maf$Hugo_Symbol, maf$Chromosome, maf$Start_Position, maf$End_Position, maf$Reference_Allele, maf$Tumor_Seq_Allele2, sep='_')
      idx.pon = which(flag.maf %in% flag.pon)
      maf = maf[-idx.pon, ]

      # Update the result table.
      sam = sub('.maf', '', maf.file)
      pth.id = unlist(strsplit(sam, '-'))[1]
      maf = cbind(batch, pth.id, sam, maf)
      names(maf)[1:3] = c('Batch', 'PTH_ID', 'Sample')
      apple = rbind(apple, maf)

      }
    }
  }

# Save the results.
write.table(apple, 'AllBatches/Result/SNV_InDel/Filtered/Variants_VariantClass_TechErrpr_PoN.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
