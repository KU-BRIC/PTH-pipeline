####################################################################################################################
####################################################################################################################
# Update referene files for filtering SNV_InDel.
# Author: Haiying Kong
# Last Modified: 13 July 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
# Get batch name and number of thread from argument passing.
while getopts ":f:l:h:t:p:n:s:c:" opt
do
  case $opt in
    r) scheme_dir="$OPTARG";;
    f) scheme_name="$OPTARG";;
    l) thresh_dp_low="$OPTARG";;
    h) thresh_dp_high="$OPTARG";;
    t) thresh_n_alt="$OPTARG";;
    p) thresh_maf_db="$OPTARG";;
    n) thresh_maf_norm="$OPTARG";;
    s) thresh_silhouette="$OPTARG";;
    c) thresh_lower_cluster_center="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done

####################################################################################################################
# Default values.
if [ -z "${scheme_dir}" ]
then
  scheme_dir=PTH
  echo "By default, the name of the filtering scheme is PTH."
fi

if [ -z "${scheme_name}" ]
then
  scheme_name=NewScheme
  echo "By default, the name of the filtering scheme is NewScheme."
fi

if [ -z "${thresh_dp_low}" ]
then 
  thresh_dp_low=200
  echo "By default, the threshold for min DP is 200."
fi

if [ -z "${thresh_dp_high}" ]
then
  thresh_dp_high=0.9975
  echo "By default, the quantile threshold for max DP is 99.75%."
fi

if [ -z "${thresh_n_alt}" ]
then 
  thresh_n_alt=4
  echo "By default, the threshold for min alternative reads is 4."
fi

if [ -z "${thresh_maf_db}" ]
then 
  thresh_maf_db=0.01
  echo "By default, to decide polymorphism, the threshold for population AF from public database is 0.01."
fi

if [ -z "${thresh_maf_norm}" ]
then 
  thresh_maf_norm=0.05
  echo "By default, to decide polymorphism, the threshold for population AF estimated from our normal samples is 0.05."
fi

if [ -z "${thresh_silhouette}" ]
then 
  thresh_silhouette=0.8
  echo "By default, to decide region specific technical error, the threshold for silhouette to decide separation of 2 clusters is 0.8."
fi

if [ -z "${thresh_lower_cluster_center}" ]
then 
  thresh_lower_cluster_center=0.08
  echo "By default, to decide region specific technical error, the threshold for the center of lower cluster is 0.08."
fi

####################################################################################################################
# Create files for the filtering scheme
filter_scheme_dir=/home/projects/cu_10184/projects/PTH/Code/Primary/FilteringScheme/${scheme_name}
Rscript ${filter_scheme_dir}/Get_Thresholds.R ${scheme_dir} ${scheme_name} ${thresh_dp_low} ${thresh_dp_high} ${thresh_n_alt} ${thresh_maf_db} ${thresh_maf_norm} ${thresh_silhouette} ${thresh_lower_cluster_center} > ${filter_scheme_dir}/Get_Thresholds.Rout 2>&1
Rscript ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError.R ${scheme_dir} ${scheme_name} > ${filter_scheme_dir}/SNP_RegionSpecificTechnicalError.Rout 2>&1

####################################################################################################################
####################################################################################################################
