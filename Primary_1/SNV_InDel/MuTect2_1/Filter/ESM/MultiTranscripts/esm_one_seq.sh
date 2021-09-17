####################################################################################################################
####################################################################################################################
# Run ESM for one variant one sequence.
# Author: Haiying Kong
# Last Modified: 30 August 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash

####################################################################################################################
# Get variable values from argument passing.
while getopts ":d:b:s:i:j:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    b) batch="$OPTARG";;
    s) sam="$OPTARG";;
    i) i="$OPTARG";;
    j) j="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done

####################################################################################################################
# Run ESM.
seq_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/temp/${sam}_${i}_${j}_seq.txt
seq=$(cat ${seq_file})

offset_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/temp/${sam}_${i}_${j}_offset.txt
offset=$(cat ${offset_file})  

in_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/temp/${sam}_${i}_${j}_input.csv
out_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/temp/${sam}_${i}_${j}_output.csv

esm_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/temp/${sam}_${i}_esm.txt

python /home/projects/cu_10184/people/haikon/Software/ESM/variant-prediction/predict.py  \
  --model-location esm1v_t33_650M_UR90S_1 esm1v_t33_650M_UR90S_2 esm1v_t33_650M_UR90S_3 esm1v_t33_650M_UR90S_4 esm1v_t33_650M_UR90S_5  \
  --sequence ${seq}  \
  --dms-input ${in_file}  \
  --mutation-col Protein_Change  \
  --dms-output ${out_file}  \
  --offset-idx ${offset}  \
  --scoring-strategy wt-marginals

cut -d',' -f 3- ${out_file} | tail -n 1 >> ${esm_file}
rm /home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/MultiTranscripts/temp/${sam}_${i}_${j}_*

####################################################################################################################
####################################################################################################################
