####################################################################################################################
####################################################################################################################
# Get list of insertions in chr13:28033736-28034557 identified by VarDict and clean them.
# Author: Haiying Kong
# Last Modified: 24 June 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
batch = args[1]
sam = args[2]
vardict.dir = args[3]
lock.dir = args[4]

####################################################################################################################
####################################################################################################################
# Read in maf column names.
maf.cols = c('Chromosome', 'Start_Position', 'End_Position', 't_alt_count', 'DP', 'AF',
             'NCBI_Build', 'Variant_Classification', 'TYPE', 'Genome_Change', 'cDNA_Change',
             'Refseq_mRNA_Id', 'Transcript_Position', 'Transcript_Strand', 'Transcript_Exon', 'Annotation_Transcript', 'Other_Transcripts',
             'Refseq_prot_Id', 'Codon_Change', 'Protein_Change',
             'Tumor_Seq_Allele2', 'SVLEN')

# Read in the list of variants from VarDict call and filter.
maf.file = paste0(vardict.dir, '/maf/', sam, '.maf')
maf = read.table(maf.file, header=TRUE, quote='', fill=TRUE, sep='\t')[ ,maf.cols]
maf = maf[(maf$Chromosome=='chr13' & maf$Start_Position>28033736 & maf$End_Position<28034557), ]
#maf = maf[-(which(maf$Variant_Type=='SNP' | maf$Variant_Type=='DEL')), ]
maf = maf[(which(maf$TYPE=='DUP' | maf$TYPE=='Insertion')), ]

# If any variants left after filtering, clean them and save.
if (nrow(maf)>0)  {
  maf$Batch = batch
  maf$Sample = sam
  maf$Length = 0
  idx = which(maf$TYPE == 'Insertion')
  if (length(idx)>0)  maf$Length[idx] = nchar(maf$Tumor_Seq_Allele2[idx])
  idx = which(maf$TYPE == 'DUP')
  if (length(idx)>0)  maf$Length[idx] = maf$SVLEN[idx]
  names(maf)[1:4] = c('Chrom', 'Start', 'End', 'N_Alt')
  maf$End = maf$Start + maf$Length - 1
  maf = maf[ ,c('Batch', 'Sample', 'Chrom', 'Start', 'End', 'Length', 'N_Alt', 'DP', 'AF', names(maf)[7:20])]
  maf = maf[order(maf$Batch,maf$Sample,maf$Chrom,maf$Start,-maf$N_Alt), ]

  # Save the result.
  write.table(maf, paste0(lock.dir, '/', sam, '.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }

####################################################################################################################
####################################################################################################################
