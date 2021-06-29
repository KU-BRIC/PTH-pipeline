####################################################################################################################
####################################################################################################################
# Update referene files for filtering SNV_InDel.
# Author: Haiying Kong
# Last Modified: 27 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

filter_scheme_dir=/home/projects/cu_10184/projects/PTH/Code/Primary/FilteringScheme/ThreeSchemes
Rscript ${filter_scheme_dir}/Get_Thresholds.R > ${filter_scheme_dir}/Get_Thresholds.Rout 2>&1
Rscript ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError.R "Long" > ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError_Long.Rout 2>&1
Rscript ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError.R "Medium" > ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError_Medium.Rout 2>&1
Rscript ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError.R "Short" > ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError_Short.Rout 2>&1

####################################################################################################################
####################################################################################################################
