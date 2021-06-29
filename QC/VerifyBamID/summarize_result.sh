#!/bin/bash

resultfile="/home/projects/cu_10184/projects/PTH/QC/Lock/VerifyBamID/Summary.txt"
echo "#SAMPLE_COMPARISON SEQ_ID RG CHIP_ID #SNPS #READS AVG_DP FREEMIX FREELK1 FREELK0 FREE_RH FREE_RA CHIPMIX CHIPLK1 CHIPLK0 CHIP_RH CHIP_RA DPREF RDPHET RDPALT" | tr ' ' '\t' > ${resultfile}

for file in /home/projects/cu_10184/projects/PTH/QC/Lock/VerifyBamID/Matched/*selfSM
do

line=$(cat ${file} | grep -v "#")
echo -n $(basename ${file}) >> ${resultfile}
echo $line | tr ' ' '\t' >> ${resultfile}

done
