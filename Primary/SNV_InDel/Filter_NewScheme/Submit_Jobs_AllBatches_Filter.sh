####################################################################################################################
####################################################################################################################
# Filter SNV-InDels.
# Author: Haiying Kong
# Last Modified: 27 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Batches:
batches=($(seq -f "%03g" 1 13))
batches=("${batches[@]/#/Primary_}")

# Submit jobs for all batches.
for batch in ${batches[@]}
do
  sh /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Filter_NewScheme/Filter.sh -d PTH -b $batch -f NewScheme
done

####################################################################################################################
####################################################################################################################
