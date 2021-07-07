 
# Automatic Pipeline for Panel Sequence Data


## Objectives.

The pipeline is developed with aims to conveniently and efficiently pre-process and analyze sequence data from targeted panels to identify somatic variants, such as SNVs, short InDels, CNVs, and internal tandem duplicates (ITDs). The pipeline also includes quality control for the sequences and produces summary statistics and plots for quality evaluations.

## Methods.

The pipeline applies existing tools on the sequence data, and works on Linux environments with cluster computer facilities.

####620110 Pre-processing.

Pre-processing of the raw sequence data applies BWA, Picard, and GATK tools to align reads on the reference genome, sort the reads, mark duplicates, and recalibrate base quality scores.

### SNV and short InDel.

SNVs and short InDels are decided as a union of variants called by VarDict, SNVer, LoFreq, and are annotated with Funcotator. SNV-InDels are filtered to exclude germline variants and variants with low pathogenic impacts.

### CNV.

CNACS and CNVkit generate two sets of CNV list, each tool excludes possible germline CNVs with panel of normals created from normal samples.

### ITD.

VarDict, Pindel, ScanITD, and getITD are applied on the reads from the targeted regions on FLT3 gene.

### QC.

FASTQuick gives summary statistics and illustrations for sequence quality evaluations.

## File organizations.

All work is saved under the folder:

    /home/projects/cu_10184/projects/PTH

### Sequence data (fastq).

The bait and target files of the panels are saved under:

    /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target
The sequence files are saved by batch:

    /home/projects/cu_10184/projects/PTH/PanelSeqData/[batch_name]/fastq

### Meta data.

Panel information for batches, and patient and tissue information for samples are saved under:

    /home/projects/cu_10184/projects/PTH/Meta
    BatchInfo.txt  SampleInfo.txt

### Reference files.

Reference files that are more specific to our study and our data, such as known variants identified from clinical study (REDCap and Horizon), filtering schemes decided by our data, panel of normals computed from our normal samples, annotation for specific regions on FLT3 gene for ITD identification, are saved under:

    /home/projects/cu_10184/projects/PTH/Reference

### Pipeline scripts.

Pipeline scripts are saved under the folder:

    /home/projects/cu_10184/projects/PTH/Code

### Files from running by a batch.

All files generated from running the pipeline on a batch are saved under:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]
Log files and error files are under:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/log
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/error
Some intermediate files and final result files are under:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Result
Quality control files are under:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/QC

## Commands and outputs.


### Full pipeline.

To run the full pipeline:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/Ensemble/Ensemble.sh -d PTH -b [batch_name] -p [panel_name] -t 8
    -d: The name of a project that is located at:
    /home/projects/cu_10184/projects/[project_name]
    -b: The name of the batch.
    -p: The name of the panel that is used for the batch. If the input is empty, it will try to find bait and target information from:
    /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target
    -t: The number of cores used by each job.

To run the full pipeline for multiple batches:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/Ensemble: vim SubmitJobs_AllBatches.sh

The outputs:
BAM files from alignment:
    
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/BAM
SNV_InDel:
Output from VarDict:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/SNV_InDel/VarDict/vcf
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/SNV_InDel/VarDict/maf
Output from SNVer:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/SNV_InDel/SNVer/vcf
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/SNV_InDel/SNVer/maf
Output from LoFreq:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/SNV_InDel/LoFreq/vcf
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/SNV_InDel/LoFreq/maf
Combined outputs from VarDict, SNVer, and LoFreq:
(Long format table has one row with one variant called by any of three callers, and can have multiple rows for the same variant if it is called by multiple callers. Wide format table reshaped Long format table, and if a variant is called by multiple callers, it takes only one row and additional information from more callers are saved in additional columns.)

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Results/SNV_InDel/AllVariants/CalledVariants_Long
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Results/SNV_InDel/AllVariants/CalledVariants_Wide
Filtered variant list:
(The combined variant list is filtered with three filtering schemes to include only potential pathogenic variants.)

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Results/SNV_InDel/Filtered/Long
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Results/SNV_InDel/Filtered/Medium
    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Results/SNV_InDel/Filtered/Short

CNV:
Output from CNVkit:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/CNV/CNVkit
Output from CNACS:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/CNV/CNACS

ITD:
Output from VarDict:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/ITD/VarDict
Output from Pindel:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/ITD/Pindel
Output from ScanITD:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/ITD/ScanITD
Output from getITD:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/ITD/getITD
IGV on FLT3 gene:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock/ITD/IGV
Combined outputs from VarDict, Pindel, ScanITD, getITD:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Result/ITD/Table
Zoom in IGV in the region of ITD:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Result/ITD/IGV
QC:
Quality evaluation on the sequence data:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/QC/FASTQuick

### Pipeline for filtering of SNV-InDel.

The main pipeline performs filtering of SNV-InDel to identify potential pathogenic variants with three fixed schemes - Long, Medium, Short.
In order to perform filtering with a new scheme with a new set of thresholds for filtering, please first create a new set of reference files for the new filtering scheme by running:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/FilteringScheme/NewScheme/Create_FilteringReferences_NewScheme.sh -r [filtering_scheme_directory] -f [filtering_scheme_name] -l [thresh_dp_low] -h [thresh_dp_high] -a [thresh_n_alt]  [thresh_maf_db] [thresh_maf_norm]  [thresh_silhouette]  [thresh_lower_cluster_center]
    -r: The name of the directory of filtering scheme. (default: PTH)
    -f: The name of the new filtering scheme. The reference files for the filtering scheme will be saved under a folder with the filtering scheme name as folder name. (default: NewScheme)
    -l: The threshold for DP to decide variants as technical error due to low DP. (default: 200)
    -h: The quantile threshold for DP to decide variants as technical error due to high DP. (default: 0.9975)
    -t: The threshold for the minimum number of alternative reads. (default: 4)
    -p: The threshold for population allele frequency obtained from public databases to decide polymorphism. (default: 0.01)
    -n: The threshold for population allele frequency estimated from our normal samples to decide polymorphism. (default: 0.05)
    -s: The threshold for silhouette score to decide if the two clusters can be considered as clearly separated when we decide region specific technical errors. (default: 0.8)
    -c: The threshold for the center of the lower cluster to decide if the lower cluster is technical error when we decide region specific technical errors. (default: 0.08)

The reference files for the new filtering scheme are saved under:

    /home/projects/cu_10184/projects/PTH/Reference/Filtering/[NewScheme]

Then, perform filtering with the new scheme by running:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Filter_NewScheme/Filter.sh -d [project_name] -b [batch_name] -f [FileringScheme]
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    -r: The name of the directory of filtering scheme. (default: PTH)
    -f: The name of the filtering scheme with reference files already created.

The variants filtered by the new scheme are saved under:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Results/SNV_InDel/Filtered/[NewScheme]

To perform filtering for multiple batches:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Filter_NewScheme/Submit_Jobs_AllBatches_Filter.sh

### CNV.

To rerun only CNVkit:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/CNV/CNVkit.sh -d [project_name] -b [batch_name]
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001

To run CNVkit for multiple batches

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/CNV/Submit_Jobs_AllBatches_CNVkit.sh

To rerun only CNACS:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/CNV/CNACS.sh -d [project_name] -b [batch_name]
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001

### ITD.

To rerun only ITD identification:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/ITD/Ensemble/ITD.sh -d [project_name] -b [batch_name] -p [panel_name] -t 8
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    -p: The name of the panel that is used for the batch. If the input is empty, it will try to find bait and target information from:
    /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target
    -t: The number of cores used by each job.
    
### QC.

To run FASTQuick:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/ITD/QC/FASTQuick.sh -d [project_name] -b [batch_name] -p [panel_name] -t 8
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    -p: The name of the panel that is used for the batch. If the input is empty, it will try to find bait and target information from:
    /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target
    -t: The number of cores used by each job.

Or to run FASTQuick for multiple batches:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/ITD/QC/SubmitJobs_AllBatches.sh

To rerun verifyBamID:

    sh /home/projects/cu_10184/projects/PTH/Code/QC/VerifyBamID/VerifyBamID.sh -d [project_name] -b [batch_name]
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    
To rerun Somalier:

    sh /home/projects/cu_10184/projects/PTH/Code/QC/Somalier/Somalier.sh -d [project_name] -b [batch_name]
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
## 
