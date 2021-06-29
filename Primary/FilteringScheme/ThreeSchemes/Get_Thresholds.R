####################################################################################################################
####################################################################################################################
# Get thresholds for 3 schemes.
# Author: Haiying Kong
# Last Modified: 27 June 2021
####################################################################################################################
####################################################################################################################
setwd('/home/projects/cu_10184/projects/PTH/')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Set values.
thresh_dp_high = 0.9975
ref.dir = '/home/projects/cu_10184/projects/PTH/Reference/Filtering'

####################################################################################################################
####################################################################################################################
# Get batch names.
batch.dir = '/home/projects/cu_10184/projects/PTH/BatchWork/'
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

####################################################################################################################
####################################################################################################################
# Long scheme:
####################################################################################################################
# Create a clean folder to save the reference files for the scheme.
filter.style = 'Long'
dir.name = paste0(ref.dir, '/', filter.style)
if (dir.exists(dir.name))  unlink(dir.name, recursive=TRUE)
dir.create(dir.name)

# Create a table for thresholds.
thresh = data.frame(scheme_name = filter.style,
                    thresh_dp_low = 200,
                    thresh_dp_high = thresh_dp,
                    thresh_dp_high_quantile = thresh_dp_high,
                    thresh_n_alt = 4,
                    thresh_maf_db = 0.02,
                    thresh_maf_norm = 0.02,
                    thresh_silhouette = 0.9,
                    thresh_lower_cluster_center = 0.05)

write.table(thresh, paste0(dir.name, '/Thresholds.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
# Medium scheme:
####################################################################################################################
# Create a clean folder to save the reference files for the scheme.
filter.style = 'Medium'
dir.name = paste0(ref.dir, '/', filter.style)
if (dir.exists(dir.name))  unlink(dir.name, recursive=TRUE)
dir.create(dir.name)

# Create a table for thresholds.
thresh = data.frame(scheme_name = filter.style,
                    thresh_dp_low = 200,
                    thresh_dp_high = thresh_dp,
                    thresh_dp_high_quantile = thresh_dp_high,
                    thresh_n_alt = 4,
                    thresh_maf_db = 0.01,
                    thresh_maf_norm = 0.05,
                    thresh_silhouette = 0.8,
                    thresh_lower_cluster_center = 0.08)

write.table(thresh, paste0(dir.name, '/Thresholds.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
# Short scheme:
####################################################################################################################
# Create a clean folder to save the reference files for the scheme.
filter.style = 'Short'
dir.name = paste0(ref.dir, '/', filter.style)
if (dir.exists(dir.name))  unlink(dir.name, recursive=TRUE)
dir.create(dir.name)

# Create a table for thresholds.
thresh = data.frame(scheme_name = filter.style,
                    thresh_dp_low = 200,
                    thresh_dp_high = thresh_dp,
                    thresh_dp_high_quantile = thresh_dp_high,
                    thresh_n_alt = 4,
                    thresh_maf_db = 0.01,
                    thresh_maf_norm = 0.05,
                    thresh_silhouette = 0.8,
                    thresh_lower_cluster_center = 0.08)

write.table(thresh, paste0(dir.name, '/Thresholds.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

####################################################################################################################
####################################################################################################################
