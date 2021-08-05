####################################################################################################################
####################################################################################################################
# Fill missing DPs and AFs in maf files from Funcotator with information from vcf.
# Author: Haiying Kong
# Last Modified: 23 July 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

# Set options and clean the space.
options(stringsAsFactors=FALSE)
rm(list=ls())

# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)
maf.dir = args[1]
sam = args[2]

####################################################################################################################
# Fill missing DP and AF in Funcotator output.
####################################################################################################################
# Read in maf file.
maf = read.table(paste0(maf.dir, '/', sam, '.maf'), header=TRUE, quote='', fill=TRUE, sep='\t')

# Fill maf with DP or AF.
maf$DP = maf$t_alt_count + maf$t_ref_count
maf$AF = maf$t_alt_count / maf$DP

####################################################################################################################
# Fix variant classification for IDH2.
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

maf = IDH2_Update(maf)

####################################################################################################################
# Update the maf file.
write.table(maf, paste0(maf.dir, '/', sam, '.maf'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
