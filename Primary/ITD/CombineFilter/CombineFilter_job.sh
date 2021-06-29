####################################################################################################################
####################################################################################################################
# Combine and filter ITDs identified from all tools.
# Author: Haiying Kong
# Last Modified: 8 June 2021
####################################################################################################################
####################################################################################################################
#!/bin/bash -i

####################################################################################################################
####################################################################################################################
# Set parameters.
scheme=Scheme_2
thresh_n_alt=1

# Get batch names.
batches=($(seq -f "%03g" 1 13))
batches=("${batches[@]/#/Primary_}")

# Submit jobs for all batches.
for batch in ${batches[@]}
do
  # Get directory names.
  Lock_ITD_dir=/home/projects/cu_10184/projects/PTH/BatchWork/${batch}/Lock/ITD
  Result_ITD_dir=/home/projects/cu_10184/projects/PTH/BatchWork/${batch}/Result/ITD
  rm -rf ${Result_ITD_dir}
  mkdir -p ${Result_ITD_dir}/Table
  mkdir ${Result_ITD_dir}/IGV

  # Get sample list.
  BAM_dir=/home/projects/cu_10184/projects/PTH/BatchWork/${batch}/Lock/BAM
  cd ${BAM_dir}
  bam_files=($(ls *.bam))
  samples=($(echo ${bam_files[@]%.bam} | tr ' ' '\n' | sort -u | tr '\n' ' '))

  for sample in ${samples[@]}
  do
    # Combine the outcome from different tools to come up with final ITD list.
    Rscript /home/projects/cu_10184/projects/PTH/Code/Source/ITD/${scheme}/Combine_Filter.R ${batch} ${sample} ${Lock_ITD_dir} ${Result_ITD_dir} ${thresh_n_alt}

    # Plot zoom in IGV if any ITD is identified for this sample.
    bat_files=${Result_ITD_dir}/IGV/${sample}*.bat

    shopt -s nullglob
    for batfile in ${Result_ITD_dir}/IGV/${sample}*.bat
    do
      rm -rf ${BAM_dir}/${sample}.bam.bai
      ln -s ${BAM_dir}/${sample}.bai ${BAM_dir}/${sample}.bam.bai

      /home/projects/cu_10184/people/haikon/Software/IGV_Linux_2.9.0/igv_auto.sh -b ${batfile}
      cat $batfile
      rm ${batfile}
      rm ${BAM_dir}/${sample}.bam.bai
    done
  done
done


####################################################################################################################
####################################################################################################################
