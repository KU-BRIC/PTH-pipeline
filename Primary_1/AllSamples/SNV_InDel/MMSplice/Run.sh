####################################################################################################################
####################################################################################################################
# Run MMSplice models for one sample.
# Author: Haiying Kong
# Last Modified: 15 November 2019
####################################################################################################################
####################################################################################################################
import numpy as np
import sys
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
# from mmsplice import predict_save
from mmsplice.utils import max_varEff

####################################################################################################################
# Set paramenters.
FastaFile = '/icgc/dkfzlsdf/analysis/B240/kong/Reference/GENCODE/GRCh37.primary_assembly.genome_nochr.fa'
GTF = '/icgc/dkfzlsdf/analysis/B240/kong/Reference/GENCODE/Homo_sapiens.GRCh37.75.chr.uniq_exon_nochr.gtf'

####################################################################################################################
####################################################################################################################
# Get arguments passed from main code.
args_values = np.delete([sys.argv], [0])
args_keys = ['tumor', 'data_dir', 'lock_dir']
args_values = args_values.tolist()
args = dict(zip(args_keys, args_values))

# Get vcf file name with full path.
vcf_file = args['data_dir'] + args['tumor'] + '.vcf'

# Define dataloader.
dl = SplicingVCFDataloader(GTF, FastaFile, vcf_file, encode=False, split_seq=True)
# dl = SplicingVCFDataloader(GTF, FastaFile, vcf_file)

# Run the MMSplice models on the tumor.
model = MMSplice()
pred = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True)
pred_max = max_varEff(pred)

# Save the results.
pred_max.to_csv(args['lock_dir'] + args['tumor'] + '.txt', index=False, sep='\t')

####################################################################################################################
####################################################################################################################
conda activate kipoi-shared__envs__kipoi-py3-keras2
