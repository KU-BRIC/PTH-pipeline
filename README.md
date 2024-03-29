 
# Automatic Pipeline for Panel Sequence Data --- Primary patient sample


## Objectives.

The pipeline is developed with aims to conveniently and efficiently pre-process and analyze sequence data from targeted panels to identify somatic variants, such as SNVs, short InDels, CNVs, and internal tandem duplicates (ITDs). The pipeline also includes quality control for the sequences and produces summary statistics and plots for quality evaluations.

## Methods.

The pipeline applies existing tools on the sequence data, and works on Linux environments with cluster computer facilities.

#### Pre-processing.

Pre-processing of the raw sequence data applies BWA, Picard, and GATK tools to align reads on the reference genome, sort the reads, mark duplicates, and recalibrate base quality scores.

#### SNV and short InDel.

SNVs and short InDels are decided as a union of variants called by VarDict, SNVer, LoFreq, and are annotated with Funcotator. SNV-InDels are filtered to exclude germline variants and variants with low pathogenic impacts.

#### CNV.

CNACS and CNVkit generate two sets of CNV list, each tool excludes possible germline CNVs with panel of normals created from normal samples.

#### ITD.

VarDict, Pindel, ScanITD, and getITD are applied on the reads from the targeted regions on FLT3 gene.

#### QC.

FASTQuick gives summary statistics and illustrations for sequence quality evaluations.

## File organizations.

All work is saved under the folder:

    /home/projects/cu_10184/projects/PTH

#### Sequence data (fastq).

The bait and target files of the panels are saved under:

    /home/projects/cu_10184/projects/[project_name]/PanelSeqData/Bait_Target

The sequence files are saved by batch:

    /home/projects/cu_10184/projects/[project_name]/PanelSeqData/[batch_name]/fastq

#### Meta data.

Panel information for batches, and patient and tissue information for samples are saved under:

    /home/projects/cu_10184/projects/[project_name]/Meta
    BatchInfo.txt  SampleInfo.txt

#### Reference files.

Reference files that are more specific to our study and our data, such as known variants identified from clinical study (REDCap and Horizon), filtering schemes decided by our data, panel of normals computed from our normal samples, annotation for specific regions on FLT3 gene for ITD identification, are saved under:

    /home/projects/cu_10184/projects/[project_name]/Reference

#### Pipeline scripts.

Pipeline scripts are saved under the folder:

    /home/projects/cu_10184/projects/PTH/Code

#### Files from running by a batch.

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


#### Full pipeline.

To run the full pipeline:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/Ensemble/Ensemble.sh -d [project_name] -b [batch_name] -p [panel_name] -t 8
    -d: The project name, for example, PTH. /home/projects/cu_10184/projects/[project_name]
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

#### SNV-InDel.
##### Call, annotation and filtering.
To call SNV_InDels with VarDict, SNVer, and LoFreq, annotate with Funcotator, and filter the variants:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Call_Anno_Filter/Call_Anno_Filter.sh -d [project_name] -b [batch_name] -p [panel_name] -t 8
    -d: The project name, for example, PTH. /home/projects/cu_10184/projects/[project_name]
    -b: The name of the batch.
    -p: The name of the panel that is used for the batch. If the input is empty, it will try to find bait and target information from:
        /home/projects/cu_10184/projects/[project_name]/PanelSeqData/Bait_Target
    -t: The number of cores used by each job.
    
To call SNV_InDels with LoFreq and annotate with Funcotator.

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/LoFreq/LoFreq.sh -d [project_name] -b [batch_name] -p [panel_name] -t 8
    -d: The project name, for example, PTH. /home/projects/cu_10184/projects/[project_name]
    -b: The name of the batch.
    -p: The name of the panel that is used for the batch. If the input is empty, it will try to find bait and target information from:
        /home/projects/cu_10184/projects/[project_name]/PanelSeqData/Bait_Target
    -t: The number of cores used by each job.

##### Filtering.
Filtering is performed to exclude technical errors and variants that are least likely pathogenically effective, such as polymorphisms or variants in intergenic regions. The process includes the following steps ("too high" or "too low" is subject to a choice of threshold):

(1) Technical errors: exclude variants identified with too low AF, or too low or too high DP.

(2) Variants classes of our interest for their high likelihood as being pathogenic: include only variants classified as the following type by Funcotator annotator.

    DE_NOVO_START_IN_FRAME, DE_NOVO_START_OUT_FRAME, Frame_Shift_Del, Frame_Shift_Ins,
    In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation,
    Splice_Site, START_CODON_SNP, Translation_Start_Site

(3) Polymorphisms: polymorphisms are decided by 

(i) variants with high population AF according to public data bases:

    dbSNP, ExAC, ESP, ClinVar, 1000G
    
(ii) variants with high population AF that is estimated by our normal samples.

(4) Pathogenic variants identified from previous studies: variants that are listed in the following databases as pathogenic are included in the final results and this overules the step (3).

    ClinVar, CIViC

(5) Region specific technical errors: if at a genetic location, more than 5 samples carry a variant, and AFs from samples form distictive 2 clusters, and the center of the lower cluster is too small, then the variants in the lower cluster are considered as technical error in a genomic region that is prone to have technical errors.

The main pipeline performs filtering of SNV-InDel to identify potential pathogenic variants with three fixed schemes - Long, Medium, Short.
The thresholds for the three schemes are saved as:

    /home/projects/cu_10184/projects/PTH/Reference/Filtering/Long/Thresholds.txt
    /home/projects/cu_10184/projects/PTH/Reference/Filtering/Medium/Thresholds.txt
    /home/projects/cu_10184/projects/PTH/Reference/Filtering/Short/Thresholds.txt

(i) For identification of polymorphisms from our normal samples, the threshold for the number of normal people who carry a given variant to decide the variant is polymorphism is decided by two different algorithms. For a given variant, the parameter estimate for population AF has point estimate and 95% confidence interval estimate. The 'Long' scheme decides the threshold under the philosophy that "Even the worst case can meet the criteria", and make the lower bound of 95% confidence interval of the parameter estimate greater than the threshold passed by the argument "-n [thresh_maf_norm]". The 'Medium' and 'Short' schemes decide the threshold under the philosophy that "The case with the highest probability can meet the criteria", and make the point estimate of population AF which has the highest probability greater than the threshold passed by the argument "-n [thresh_maf_norm]".

(ii) For the list for exclusion candidate, the 'Long' and 'Medium' schemes exclude variants that do not have any information on the column "COSMIC_tissue_types_affected" from Funcotator annotation, while the 'Short' scheme excludes variants that do not have any information on the column "COSMIC_overlapping_mutations" from Funcotator annotation.

(6) White list and black list.
white list are saved as:

    /home/projects/cu_10184/projects/PTH/Reference/Filtering/BlackList/white.bed
    /home/projects/cu_10184/projects/PTH/Reference/Filtering/BlackList/white.vcf

black list are saved as:

    /home/projects/cu_10184/projects/PTH/Reference/Filtering/BlackList/black.bed
    /home/projects/cu_10184/projects/PTH/Reference/Filtering/BlackList/black.vcf

If there is contradiction between white and black list, black list wins.

###### Three filtering schemes.

In order to perform filtering with three schemes that are installed in the pipeline, please run:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/SNV_InDel/Filter_ThreeSchemes/Filter.sh -d [project_name] -b [batch_name]

###### New filtering scheme.

In order to perform filtering with a new scheme that is modified from the 'Medium' scheme with a new set of thresholds, please first create a new set of reference files for the new filtering scheme by running:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/FilteringScheme/NewScheme/Create_FilteringReferences_NewScheme.sh -r [filtering_scheme_directory] -f [filtering_scheme_name] -l [thresh_dp_low] -h [thresh_dp_high] -t [thresh_n_alt] -p [thresh_maf_db] -n [thresh_maf_norm]  -s [thresh_silhouette] -c  [thresh_lower_cluster_center]
    -r: The name of the directory of filtering scheme. (default: PTH)
    -f: The name of the new filtering scheme. The reference files for the filtering scheme will be saved under a folder with the filtering scheme name as folder name. (default: NewScheme)
    -l: The threshold for DP to decide variants as technical error due to low DP. (default: 200)
    -h: The quantile threshold for DP to decide variants as technical error due to high DP. (default: 0.9975)
    -t: The threshold for the minimum number of alternative reads. (default: 4)
    -p: The threshold for population allele frequency obtained from public databases to decide polymorphism. (default: 0.01)
    -n: The threshold for population allele frequency estimated from our normal samples to decide polymorphism. (default: 0.05)
    -s: The threshold for silhouette score to decide if the two clusters can be considered as clearly separated when we decide region specific technical errors. (default: 0.8)
    -c: The threshold for the cente/home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Targetr of the lower cluster to decide if the lower cluster is technical error when we decide region specific technical errors. (default: 0.08)

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

#### CNV.

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

#### ITD.

To rerun only ITD identification:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/ITD/Ensemble/ITD.sh -d [project_name] -b [batch_name] -p [panel_name] -t 8
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    -p: The name of the panel that is used for the batch. If the input is empty, it will try to find bait and target information from:
    /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target
    -t: The number of cores used by each job.
    
#### QC.

To run FASTQuick:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/QC/FASTQuick.sh -d [project_name] -b [batch_name] -p [panel_name] -t 8
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    -p: The name of the panel that is used for the batch. If the input is empty, it will try to find bait and target information from:
    /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target
    -t: The number of cores used by each job.

Or to run FASTQuick for multiple batches:

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/QC/SubmitJobs_AllBatches.sh

To rerun verifyBamID:

    sh /home/projects/cu_10184/projects/PTH/Code/QC/VerifyBamID/VerifyBamID.sh -d [project_name] -b [batch_name]
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    
To rerun Somalier:

    sh /home/projects/cu_10184/projects/PTH/Code/QC/Somalier/Somalier.sh -d [project_name] -b [batch_name]
    -d: The project name, for example, PTH.
    -b: The batch name, for example, Primary_001
    
To run only the part of pipeline to get our QC summary table:

    sh /home/projects/cu_10184/projects/PTH/Code/QC/QC_Summary/SubmitJobs_AllBatches.sh

Then, the QC summary table of each sample can be found as:

    /home/projects/cu_10184/projects/PTH/QC/Result/FASTQuick/ByBatch/[batch_name]/[sample_name]/OurSummary.txt

To get QC summary table for all samples from all batches:

    Rscript /home/projects/cu_10184/projects/PTH/Code/QC/QC_Summary/AllSamples_QC_Summary.R

The result will be at:

    /home/projects/cu_10184/projects/PTH/QC/Result/FASTQuick/AllBatches/QC_Summary.txt

#### To run pipeline with a new panel.
(1) Save target files as bed files under the designated directories.

   The target file before padding, with chromosomes denoted as: chr1, chr2, ..., chr22, chrX, chrY:
    
    /home/projects/cu_10184/projects/[project_name]/PanelSeqData/Bait_Target/Chr_Original

   The target file with regions padded, with chromosomes denoted as: chr1, chr2, ..., chr22, chrX, chrY:
    
    /home/projects/cu_10184/projects/[project_name]/PanelSeqData/Bait_Target/Chr_Padded
   
   The target file with regions padded, with chromosomes denoted as: 1, 2, ..., 22, X, Y:
    
    /home/projects/cu_10184/projects/[project_name]/PanelSeqData/Bait_Target/Padded

(2) Create or update Meta file.

    /home/projects/cu_10184/projects/[project_name]/Meta/BatchInfo.txt
    
   The file should be tab deliminated file in the format:
   
    Batch_ID          Bait    Target
    [name_of_batch]   [na]    [target_file_name_without_path_without_.bed]
    
   Bait file is not called by pipeline, and therefore can be fill with any string for the cell.
   
(3) Run the pipeline without passing any values to the argument -p. For instance,

    sh /home/projects/cu_10184/projects/PTH/Code/Primary/Ensemble/Ensemble.sh -d PTH -b [batch_name] -t 8


# Automatic Pipeline for Panel Sequence Data --- PDXample

## Classification of reads into 5 categories: human, mouse, both, neither, ambiguous.
Summary table for the reads in each category can be found in:

    /home/projects/cu_10184/projects/PTH/PanelSeqData/PDX_001/meta/Summary_ReadCounts_xengsort.txt

## Reads from human are analyzed to identify genetic variants with the command:

    sh /home/projects/cu_10184/projects/PTH/Code/PDX/Ensemble/Ensemble.sh -d PTH -b [batch_name] -t [num_thread] -c xengsort -h [thresh_num_human_reads]
    -b: batch name.
    -t: number of threads to use on computerome.
    -h: minimum number of reads from human to be analyzed.

### Intermediate results can be found:

    /home/projects/cu_10184/projects/PTH/BatchWork/[batch_name]/Lock
    
### Final results can be found:

    /home/projects/cu_10184/projects/PTH/BatchWork/PDX_001/Result
    

# Comprehensive analysis on primary and PDX samples.

## The pedigree information for primary and PDX samples are saved with the fixed format as:

    /home/projects/cu_10184/projects/PTH/Meta/PAT_PDX_Info.txt

## The target file for PDX samples is saved as:

    /home/projects/cu_10184/projects/PTH/PanelSeqData/Bait_Target/Target/Chr_Padded/Focused_myeloid_panel-All_target_segments_covered_by_probes-TE-93310852_hg38_v2_190722165759.bed

For primary samples, variants are filtered for PDX target regions.

## To run the analysis:

    Rscript /home/projects/cu_10184/projects/PTH/Code/PAT_PDX/Evolution.R [thresh_spaghetti] > /home/projects/cu_10184/projects/PTH/Code/PAT_PDX/Evolution.Rout
    [thresh_spaghetti]: On spagjetti plot, only variants that have AF greater than [thresh_spaghetti] either in PAT or PDX sample will be on the plot.
    
## The results are under the directory:

    /home/projects/cu_10184/projects/PTH/PAT_PDX/Result

##
