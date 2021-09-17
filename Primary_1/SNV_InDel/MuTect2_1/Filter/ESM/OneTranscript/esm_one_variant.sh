####################################################################################################################
####################################################################################################################
# Run ESM for one variant.
# Author: Haiying Kong
# Last Modified: 26 August 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash

####################################################################################################################
# Get variable values from argument passing.
while getopts ":d:b:s:i:" opt
do
  case $opt in
    d) dir_name="$OPTARG";;
    b) batch="$OPTARG";;
    s) sam="$OPTARG";;
    i) i="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done

####################################################################################################################
# Run ESM.
seq_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/temp/${sam}_${i}_seq.txt
seq=$(cat ${seq_file})

offset_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/temp/${sam}_${i}_offset.txt
offset=$(cat ${offset_file})  

acorn_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/temp/${sam}_${i}_acorn.maf
in_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/temp/${sam}_${i}_input.csv
out_csv_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/temp/${sam}_${i}_output.csv
out_tab_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/temp/${sam}_${i}_output.tab
res_file=/home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/Original/${sam}_SAAS.maf

python /home/projects/cu_10184/people/haikon/Software/ESM/variant-prediction/predict.py  \
  --model-location esm1v_t33_650M_UR90S_1 esm1v_t33_650M_UR90S_2 esm1v_t33_650M_UR90S_3 esm1v_t33_650M_UR90S_4 esm1v_t33_650M_UR90S_5  \
  --sequence ${seq}  \
  --dms-input ${in_file}  \
  --mutation-col Protein_Change  \
  --dms-output ${out_csv_file}  \
  --offset-idx ${offset}  \
  --scoring-strategy wt-marginals

cut -d',' -f 3- ${out_csv_file} | tail -n 1 | sed 's/\,/\t/g' > ${out_tab_file}
paste ${acorn_file} ${out_tab_file} >> ${res_file}
rm /home/projects/cu_10184/projects/${dir_name}/BatchWork_1/${batch}/Result/SNV_InDel/MuTect2_1/ESM/OneTranscript/temp/${sam}_${i}_*

####################################################################################################################
####################################################################################################################
