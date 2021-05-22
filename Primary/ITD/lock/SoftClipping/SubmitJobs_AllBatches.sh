####################################################################################################################
####################################################################################################################
# Call ITD with soft-clipping for all batches.
# Author: Haiying Kong
# Last Modified: 28 March 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Batches:
batches=(Batch001 Batch002 Batch003 Batch004 Batch005 Batch006 Batch007 Batch008 Batch009 Batch010 Batch011)

# Submit jobs for all batches.
for batch in ${batches[@]}
do
  sh /home/projects/cu_10145/people/haikon/Project/PTH/Code/ITD/SoftClipping/SoftClipping.sh -b $batch
done


####################################################################################################################
####################################################################################################################
