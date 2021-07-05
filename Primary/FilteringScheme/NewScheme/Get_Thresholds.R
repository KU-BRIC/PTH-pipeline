####################################################################################################################
####################################################################################################################
# Get threshold for the filtering scheme.
# Author: Haiying Kong
# Last Modified: 3 July 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Get argument values from command line.
args = commandArgs(trailingOnly=TRUE)

scheme_dir = args[1]
scheme_name = args[2]
thresh_dp_low = as.numeric(args[3])
thresh_dp_high = as.numeric(args[4])
thresh_n_alt = as.numeric(args[5])
thresh_maf_db = as.numeric(args[6])
thresh_maf_norm = as.numeric(args[7])
thresh_silhouette = as.numeric(args[8])
thresh_lower_cluster_center = as.numeric(args[9])

####################################################################################################################
####################################################################################################################
# Define directory to save results.
ref.dir = paste0('/home/projects/cu_10184/projects/', scheme_dir, '/Reference/Filtering/', scheme_name)
if (dir.exists(ref.dir))   unlink(ref.dir, recursive=TRUE)
dir.create(ref.dir, recursive=TRUE)

####################################################################################################################
####################################################################################################################
# Get batch names.
batch.dir = '/home/projects/cu_10184/projects/', scheme_dir, '/BatchWork/'
batches = dir(batch.dir, pattern='^Primary_')

# Collect all DPs from all samples.
DP = c()
for (batch in batches)   {
  sam.dir = paste0(batch.dir, batch, '/Result/SNV_InDel/AllVariants/Callers_Wide/')
  samples = dir(sam.dir, pattern='.maf')
  for (sam in samples)   {
    dp = read.table(paste0(sam.dir, sam), header=TRUE, colClasses=c(rep('NULL',18), 'integer', rep('NULL',98)), quote='', fill=TRUE, sep='\t')$DP
    DP = c(DP, dp)
    }
  }

# Get threshold for high DP.
thresh_dp = quantile(DP, thresh_dp_high)

# Create a table for thresholds.
thresh = data.frame(scheme_name = scheme_name,
                    thresh_dp_low = thresh_dp_low,
                    thresh_dp_high = thresh_dp,
                    thresh_dp_high_quantile = thresh_dp_high,
                    thresh_n_alt = thresh_n_alt,
                    thresh_maf_db = thresh_maf_db,
                    thresh_maf_norm = thresh_maf_norm,
                    thresh_silhouette = thresh_silhouette,
                    thresh_lower_cluster_center = thresh_lower_cluster_center)

write.table(thresh, paste0(ref.dir, '/Thresholds.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
