####################################################################################################################
####################################################################################################################
# Update referene files for filtering with all data.
# Author: Haiying Kong
# Last Modified: 27 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

filter_scheme_dir=/home/projects/cu_10184/projects/PTH/Code/Primary/FilteringScheme
R CMD BATCH ${filter_scheme_dir}/ThreeSchemes/Get_thresh_HighDP.R
R CMD BATCH ${filter_scheme_dir}/ThreeSchemes/Long_SNP_RegionSpecificTechnicalError.R
R CMD BATCH ${filter_scheme_dir}/ThreeSchemes/Medium_SNP_RegionSpecificTechnicalError.R
R CMD BATCH ${filter_scheme_dir}/ThreeSchemes/Short_SNP_RegionSpecificTechnicalError.R

####################################################################################################################
####################################################################################################################
