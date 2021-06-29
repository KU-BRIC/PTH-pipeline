####################################################################################################################
####################################################################################################################
# Get list of Horizon and REDCap variants called by our pipeline.
# Author: Haiying Kong
# Last Modified: 3 June 2021
####################################################################################################################
####################################################################################################################
#!/home/projects/cu_10184/people/haikon/Software/R-4.0.4/bin/Rscript

setwd('/home/projects/cu_10184/projects/PTH')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
####################################################################################################################
# Read in table with information for AF from our calls and reference data.
af = read.table('AllBatches/SNV_InDel/ClinicTag/CalledReferenceVariants.txt', header=TRUE, sep='\t')
af = af[ ,grep('AF', names(af))]

####################################################################################################################
# PLots.
####################################################################################################################
pdf('AllBatches/SNV_InDel/ClinicTag/Plot_AF_Called_Reference.pdf')

# AF_Clinic v.s. AF:
plot(af$AF_Clinic, af$AF, pch=20, col='dodgerblue1',
     main='AF association', xlab='AF_Clinic', ylab='AF_max')
m = lm(AF~AF_Clinic, data=af)
abline(m, col='brown2')
correlation = round(cor(af$AF_Clinic, af$AF, use='complete.obs', method='pearson'), 3)
legend(0.1, 0.85, legend=paste0('Correlation Coefficient: ', correlation), col=c('green'), cex=0.8, bty='n')

# AF_Clinic v.s. AF_VarDict:
plot(af$AF_Clinic, af$AF_VarDict, pch=20, col='dodgerblue1',
     main='AF association', xlab='AF_Clinic', ylab='AF_VarDict')
m = lm(AF_VarDict~AF_Clinic, data=af)
abline(m, col='brown2')
correlation = round(cor(af$AF_Clinic, af$AF_VarDict, use='complete.obs', method='pearson'), 3)
legend(0.1, 0.85, legend=paste0('Correlation Coefficient: ', correlation), col=c('green'), cex=0.8, bty='n')

# AF_Clinic v.s. AF_SNVer:
plot(af$AF_Clinic, af$AF_SNVer, pch=20, col='dodgerblue1',
     main='AF association', xlab='AF_Clinic', ylab='AF_SNVer')
m = lm(AF_SNVer~AF_Clinic, data=af)
abline(m, col='brown2')
correlation = round(cor(af$AF_Clinic, af$AF_SNVer, use='complete.obs', method='pearson'), 3)
legend(0.1, 0.85, legend=paste0('Correlation Coefficient: ', correlation), col=c('green'), cex=0.8, bty='n')

# AF_Clinic v.s. AF_LoFreq:
plot(af$AF_Clinic, af$AF_LoFreq, pch=20, col='dodgerblue1',
     main='AF association', xlab='AF_Clinic', ylab='AF_LoFreq')
m = lm(AF_LoFreq~AF_Clinic, data=af)
abline(m, col='brown2')
correlation = round(cor(af$AF_Clinic, af$AF_LoFreq, use='complete.obs', method='pearson'), 3)
legend(0.1, 0.85, legend=paste0('Correlation Coefficient: ', correlation), col=c('green'), cex=0.8, bty='n')

dev.off()

####################################################################################################################
####################################################################################################################
