####################################################################################################################
####################################################################################################################
# Create large table for all variants with ESM scores.
# Author: Haiying Kong
# Last Modified: 5 September 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH/temp')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(stringr)

####################################################################################################################
# Set parameters.
proj.name = 'PTH'
batches = paste0('Primary_', str_pad(1:13, 3, pad='0'))

####################################################################################################################
# Read in the names of columns that will be kept.
maf.cols = c('Batch', 'Sample', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Reference_Allele', 'Tumor_Seq_Allele2',
             'Refseq_mRNA_Id', 'Genome_Change', 'cDNA_Change', 'Protein_Change', 'Variant_Classification', 'Variant_Type',
             'COSMIC_tissue_types_affected', 'COSMIC_overlapping_mutations', 'AF', 'DP')
maf.cols.esm = c(maf.cols, paste0('esm1v_t33_650M_UR90S_', 1:5))

####################################################################################################################
####################################################################################################################
# Create large maf file by concatenating all maf files from all samples.
####################################################################################################################
# Collect all variants from all samples with technical error, variant class and PoN filtered.
SAAS = c()
rest = c()

for (batch in batches)  {
  maf.dir = paste0('/home/projects/cu_10184/projects/PTH/BatchWork_1/', batch, '/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Original')

  # Single amino acid substitution:
  maf.files = dir(maf.dir, pattern='_SAAS.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')    
    if (nrow(maf)>0)  {
      # Add and trim columns.
      maf = cbind(batch, sub('_SAAS.maf', '', maf.file), maf)
      names(maf)[1:2] = c('Batch', 'Sample')
      maf = maf[ ,maf.cols.esm]

      # Update apple.
      SAAS = rbind(SAAS, maf)
      }
    }

  # Rest:
  maf.files = dir(maf.dir, pattern='_rest.maf')
  for (maf.file in maf.files)  {
    maf = read.table(paste0(maf.dir, '/', maf.file), header=TRUE, quote='', fill=TRUE, sep='\t')
    if (nrow(maf)>0)  {
      # Add and trim columns.
      maf = cbind(batch, sub('_rest.maf', '', maf.file), maf)
      names(maf)[1:2] = c('Batch', 'Sample')
      maf = maf[ ,maf.cols]
      
      # Update apple.
      rest = rbind(rest, maf)
      }
    }
  }

# Save the results.
write.table(SAAS, '/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Apple_SAAS.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.table(rest, '/home/projects/cu_10184/projects/PTH/AllBatches_1/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/Apple_rest.txt',
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
