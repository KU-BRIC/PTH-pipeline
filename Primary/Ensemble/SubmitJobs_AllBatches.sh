####################################################################################################################
####################################################################################################################
# Pre-process, call variants, and filter.
# Author: Haiying Kong
# Last Modified: 10 April 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Batches:
batches=(Primary_001 Primary_002 Primary_003 Primary_004 Primary_005 Primary_006 Primary_007 Primary_008 Primary_009 Primary_010 Primary_011)

# Submit jobs for all batches.
for batch in ${batches[@]}
do
  sh /home/projects/cu_10184/projects/PTH/Code/Primary/Ensemble/Ensemble.sh -d PTH -b $batch -t 8
done


####################################################################################################################
####################################################################################################################