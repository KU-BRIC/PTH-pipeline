####################################################################################################################
####################################################################################################################
# Create PoN with Step2 and Steo3 for MuTect2.
# Author: Haiying Kong
# Last Modified: 19 July 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

dir_name=PTH
log_file=/home/projects/cu_10184/projects/${dir_name}/AllBatches_1/log/SNV_InDel/MuTect2/PoN_Step2_Step3.log
error_file=/home/projects/cu_10184/projects/${dir_name}/AllBatches_1/error/SNV_InDel/MuTect2/PoN_Step2_Step3.error

qsub -o ${log_file} -e ${error_file} -N MuTect2_PoN_Step2_Step3  \
  -l nodes=1:ppn=8,mem=70gb,walltime=70:00:00  \
  /home/projects/cu_10184/projects/PTH/Code/Primary_1/SNV_InDel/MuTect2/PoN/Step2_Step3_Create_GenomicsDB_PoN_job.sh

####################################################################################################################
####################################################################################################################
