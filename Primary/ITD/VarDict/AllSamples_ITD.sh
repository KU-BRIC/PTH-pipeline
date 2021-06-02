####################################################################################################################
####################################################################################################################
# Identify ITD with VarDict.
# Author: Haiying Kong
# Last Modified: 1 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Get batch names.
batches=($(seq -f "%03g" 1 13))
batches=("${batches[@]/#/Primary_}")

# Run cleaning ScanITD for all samples 
for batch in ${batches[@]}
do
  maf_dir=/home/projects/cu_10184/data/projects/PTH/BatchWork/$batch/Lock/SNV_InDel/VarDict
  itd_dir=/home/projects/cu_10184/data/projects/PTH/BatchWork/$batch/Lock/ITD/VarDict

  # Get list of samples.
  cd ${maf_dir}/maf
  maf_files=($(ls *.maf))
  samples=($(echo ${maf_files[@]%.maf} | tr ' ' '\n' | sort -u | tr '\n' ' '))

  for sample in ${samples[@]}
  do
    Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/Scheme_2/VarDict_FilterClean.R ${batch} ${sample} ${maf_dir} ${itd_dir}
  done
done

####################################################################################################################
####################################################################################################################
