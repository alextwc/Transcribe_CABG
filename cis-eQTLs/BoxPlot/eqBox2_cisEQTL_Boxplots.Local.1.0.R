#!/usr/bin/Rscript
message("############################################################################################################")
message("## R-script to run enhanced eqBox2() to automatically generate the boxplots of cis-eQTLs results          ##")
message("##                                                                                                        ##")
message("##  Author: Dr. Alex Tzuu-Wang Chang                                                                      ##")
message("##    Date: 2016-09-23                                                                                    ##")
message("##                                                                                                        ##")
message("## Version: 1.0 Setting the path variables for imputed Transcribe-118 genotyping dataset                  ##")
message("##                                                                                                        ##") 
message("## Require: R v3.3.1, python etc.                                                                         ##")
message("##                                                                                                        ##")
message("## To install Bioconductor: source(\"http://bioconductor.org/biocLite.R\")                                  ##")
message("##                          biocLite(\"BiocInstaller\")                                                     ##")
message("## To upgrade Bioconductor: source(\"http://bioconductor.org/biocLite.R\")                                  ##")
message("##                          biocLite(\"BiocUpgrade\")                                                       ##")
message("## To install a package, for example the \"SummarizedExperiment\". Run the following commands:              ##")
message("##                          source(\"http://bioconductor.org/biocLite.R\")                                  ##")
message("##                          biocLite(\"SummarizedExperiment\")                                              ##")
message("##                                                                                                        ##")
message("## When submit jobs at the outside of Channing:                                                           ##")
message("## (1) Login to capecod first (2) qsub -l lx yourScript.sh                                                ##")
message("## Switch to LINUX node from outside: ssh -l retwc -o StrictHostKeyChecking=no nantucket.bwh.harvard.edu  ##")
message("## Switch to LINUX node from  inside: ssh -l retwc nantucket                                              ##")
message("## Switch to ALKAN node from  inside: ssh alkan01                                                         ##")
message("## alias qs1='qrsh -pe ompi 8 -l virtual_free=32G -l m_core=5'                                            ##")
message("## alias qs2='qrsh -pe ompi 8 -l virtual_free=92G -l m_core=9'                                            ##")
message("##                                                                                                        ##")
message("## Run R in Channing:[ qs | (qrsh -l lx6,large) | (qrsh -l lx6,12hour=true) | (qrsh -l lx) ] then R;      ##")
message("## qrsh -l lx6 -l mem_free=92G -l m_core=13; qrsh -l lx6 -l mem_free=80G -l m_core=5;                     ##")
message("## qrsh -l rstudio; R --vanilla; /udd/stvjc/VM/R-devel-dist/bin/R; ~stvjc/bin/Rdevel --vanilla;           ##")
message("## `/udd/stvjc/bin/Rdevel RHOME`/bin/Rscript                                                              ##")
message("##                                                                                                        ##")
message("## Run R script in Channing:                                                                              ##")
message("## nohup R --no-save < myRscript.R > myRscriptR.out &                                                     ##")
message("## R --vanilla < /udd/retwc/R/test.R > /udd/retwc/R/testR.out ; qsub -l lx ./runMyRscript.sh              ##")
message("##                                                                                                        ##")
message("## After qs or qrsh the R executable file is located at the following path:                               ##")
message("## /local/bin/R -> /app/R-3.3.1@i86-rhel6.0/bin/R                                                         ##")
message("## Rscript -> /app/R-3.3.1@i86-rhel6.0/bin/Rscript                                                        ##")
message("## Check Job status: qds retwc; qstat -r; qstat -f; qstat -f | grep retwc | less; qstat -ls; qhost;       ##")
message("##                                                                                                        ##")
message("## Usage1 (for Channing Cluster): ( Rscript ./test.R ./LV135preX.csv ./LV133postX.csv ) >& testR.log      ##")
message("## Usage2 (for Channing Cluster): qbR -A \"./LV135preX.csv ./LV133postX.csv\" ./test.R                      ##")
message("## Result: Usage 1 works; Usage 2 (qbR command line) is also working now after adding -A                  ##")
message("############################################################################################################")
print(Sys.time())

R.version
R.home()
.libPaths()
sessionInfo()
# installed.packages()
# .libPaths("/udd/retwc/R/library/3.1")
# .libPaths( c( .libPaths(), "/udd/retwc/R/library/3.1") )
# detach("package:GGtools" , unload=TRUE)
# install.packages("name-of-your-package", lib="~/R/library")
# base.dir <- getwd();
#     odir <- getwd();
message("###############################################################################")
message("# Environmental Variables and Working directory settings, R packages loading  #")
message("###############################################################################")
print(Sys.time())
# require("peer", lib.loc="/udd/retwc/R")
require("VariantAnnotation")
require("GenomicRanges")
require("GenomeInfoDb")
require("data.table")
require("Rsamtools")
require("gQTLstats")
require("gQTLBase")
require("ggplot2")
require("reshape2")
require("graphics")
require("biomaRt")
require("tools")
listMarts()
set.seed(12345)
d <- Sys.Date()
sourceDir <- "D:/BWHMS/Meetings02/2016-09-22";
source(paste(sourceDir, "/InputFiles/eeqBox2.R", sep = "")); # Loading the enhanced eqBox2() function
                                                             # The enhanced eqBox2 function named by "eeqBox2"
setwd("D:/BWHMS/Meetings02/2016-09-29")
imputed_all <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_AllChr22_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_01  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr01_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_02  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr02_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_03  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr03_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_04  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr04_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_05  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr05_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_06  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr06_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_07  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr07_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_08  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr08_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_09  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr09_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_10  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr10_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_11  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr11_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_12  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr12_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_13  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr13_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_14  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr14_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_15  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr15_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_16  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr16_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_17  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr17_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_18  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr18_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_19  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr19_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_20  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr20_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_21  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr21_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_22  <- "D:/BWHMS/TRANSCRiBE118Genotypes/VCF-IR/Transcribe118_Chr22_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_23  <- "D:/BWHMS/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh23.vcf-IR.vcf.gz"
vcfpath_24  <- "D:/BWHMS/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh24.vcf-IR.vcf.gz"
vcfpath_25  <- "D:/BWHMS/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh25.vcf-IR.vcf.gz"
vcfpath_26  <- "D:/BWHMS/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh26.vcf-IR.vcf.gz"
tf01 = TabixFile(vcfpath_01)
tf02 = TabixFile(vcfpath_02)
tf03 = TabixFile(vcfpath_03)
tf04 = TabixFile(vcfpath_04)
tf05 = TabixFile(vcfpath_05)
tf06 = TabixFile(vcfpath_06)
tf07 = TabixFile(vcfpath_07)
tf08 = TabixFile(vcfpath_08)
tf09 = TabixFile(vcfpath_09)
tf10 = TabixFile(vcfpath_10)
tf11 = TabixFile(vcfpath_11)
tf12 = TabixFile(vcfpath_12)
tf13 = TabixFile(vcfpath_13)
tf14 = TabixFile(vcfpath_14)
tf15 = TabixFile(vcfpath_15)
tf16 = TabixFile(vcfpath_16)
tf17 = TabixFile(vcfpath_17)
tf18 = TabixFile(vcfpath_18)
tf19 = TabixFile(vcfpath_19)
tf20 = TabixFile(vcfpath_20)
tf21 = TabixFile(vcfpath_21)
tf22 = TabixFile(vcfpath_22)
tf23 = TabixFile(vcfpath_23)
tf24 = TabixFile(vcfpath_24)
tf25 = TabixFile(vcfpath_25)
tf26 = TabixFile(vcfpath_26)

rD_Exclude3PCs <- get(load(paste(sourceDir, "/ResponseQTL/Without3PCs/imputed_Delta/ResponseQTL_Delta_imputed_deltaResidual_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)
rD_Include3PCs <- get(load(paste(sourceDir,    "/ResponseQTL/With3PCs/imputed_Delta/ResponseQTL_Delta_imputed_deltaResidual_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)
rR_Exclude3PCs <- get(load(paste(sourceDir, "/ResponseQTL/Without3PCs/imputed_Ratio/ResponseQTL_Ratio_imputed_ratioResidual_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)
rR_Include3PCs <- get(load(paste(sourceDir,    "/ResponseQTL/With3PCs/imputed_Ratio/ResponseQTL_Ratio_imputed_ratioResidual_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)
sB_Exclude3PCs <- get(load(paste(sourceDir, "/StateseQTL/Without3PCs/imputed_Baseline/StateseQTL_Baseline_imputed_exprPre_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)
sB_Include3PCs <- get(load(paste(sourceDir,    "/StateseQTL/With3PCs/imputed_Baseline/StateseQTL_Baseline_imputed_exprPre_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)
sI_Exclude3PCs <- get(load(paste(sourceDir, "/StateseQTL/Without3PCs/imputed_Ischemia/StateseQTL_Ischemia_imputed_exprPost_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)
sI_Include3PCs <- get(load(paste(sourceDir,    "/StateseQTL/With3PCs/imputed_Ischemia/StateseQTL_Ischemia_imputed_exprPost_FDR.05_PEER-10.rda", sep = "")))
rm(exprSeQTL_fdr.05)

EXP <- list(df1 = rD_Exclude3PCs, df2 = rR_Exclude3PCs, df3 = rD_Include3PCs, df4 = rR_Include3PCs)
EXP <- lapply(EXP, function(df) {
	 df <- df[sample(1:nrow(df)), ]
	 df <- df[with(df, order(-chisq, piFDR, mindist)), ]
	 df <- df[!duplicated(df$probeid), ]
	 df <- df[order(df$chisq, decreasing = TRUE), ]
	 df$pValue = pchisq(df$chisq, 1, lower.tail=FALSE)
	 df
})
rD_Exclude3PCs <- EXP$df1
rR_Exclude3PCs <- EXP$df2
rD_Include3PCs <- EXP$df3
rR_Include3PCs <- EXP$df4
rm(EXP)

EXP <- list(df1 = sB_Exclude3PCs, df2 = sI_Exclude3PCs, df3 = sB_Include3PCs, df4 = sI_Include3PCs)
EXP <- lapply(EXP, function(df) {
	 df <- df[sample(1:nrow(df)), ]
	 df <- df[with(df, order(-chisq, piFDR, mindist)), ]
	 df <- df[!duplicated(df$probeid), ]
	 df <- df[order(df$chisq, decreasing = TRUE), ]
	 df$pValue = pchisq(df$chisq, 1, lower.tail=FALSE)
	 df
})
sB_Exclude3PCs <- EXP$df1
sI_Exclude3PCs <- EXP$df2
sB_Include3PCs <- EXP$df3
sI_Include3PCs <- EXP$df4
rm(EXP)

for (i in c("rD_Exclude3PCs","rR_Exclude3PCs","rD_Include3PCs","rR_Include3PCs")){
	pdf(file=paste("piFDR_pValue_", i, "_ResponseQTL.pdf", sep=""))
	dat <- get(i)
	n   <- nrow(dat)
	plot(dat$pValue, dat$piFDR, type = "p", col = "red", pch=20, cex=1.0, main = paste("piFDR vs. p-Value of ", i, " eQTLs", sep=""),
	xlab = paste("pValues of ", i, " most significant ", n, " SNPs", sep=""), 
	ylab = paste("piFDR of ", i, " most significant ", n, " SNPs", sep=""))
  dev.off()
}

for (i in c("rD_Exclude3PCs","rR_Exclude3PCs","rD_Include3PCs","rR_Include3PCs")){
	pdf(file=paste("piFDR_adjustedPvalue_", i, "_ResponseQTL.pdf", sep=""))
	dat <- get(i)
	n   <- nrow(dat)
	plot(p.adjust(dat$pValue, "BH"), dat$piFDR, type = "p", col = "red4", pch=20, cex=1.0, main = paste("piFDR vs. adjusted p-Value of ", i, " eQTLs", sep=""),
	xlab = paste("BH adjusted pValues of ", i, " most significant ", n, " SNPs", sep=""), 
	ylab = paste("piFDR of ", i, " most significant ", n, " SNPs", sep=""))
  dev.off()
}

for (i in c("rD_Exclude3PCs","rR_Exclude3PCs","rD_Include3PCs","rR_Include3PCs")){
	dat <- get(i)
	n   <- nrow(dat)
  write.csv(dat, file=paste("./ResponseQTL_",i, "_FDR.05_PEER10_", n, "_topSNPs_", d, ".csv", sep=""), row.names=FALSE)
}

for (i in c("sB_Exclude3PCs","sI_Exclude3PCs","sB_Include3PCs","sI_Include3PCs")){
	pdf(file=paste("piFDR_pValue_", i, "_StateseQTL.pdf", sep=""))
	dat <- get(i)
	n   <- nrow(dat)
	plot(dat$pValue, dat$piFDR, type = "p", col = "red", pch=20, cex=0.4, main = paste("piFDR vs. p-Value of ", i, " eQTLs", sep=""),
	xlab = paste("pValues of ", i, " most significant ", n, " SNPs", sep=""), 
	ylab = paste("piFDR of ", i, " most significant ", n, " SNPs", sep=""))
  dev.off()
}

for (i in c("sB_Exclude3PCs","sI_Exclude3PCs","sB_Include3PCs","sI_Include3PCs")){
	pdf(file=paste("piFDR_adjustedPvalue_", i, "_StateseQTL.pdf", sep=""))
	dat <- get(i)
	n   <- nrow(dat)
	plot(p.adjust(dat$pValue, "BH"), dat$piFDR, type = "p", col = "red4", pch=20, cex=0.4, main = paste("piFDR vs. adjusted p-Value of ", i, " eQTLs", sep=""),
	xlab = paste("BH adjusted pValues of ", i, " most significant ", n, " SNPs", sep=""),
	ylab = paste("piFDR of ", i, " most significant ", n, " SNPs", sep=""))
  dev.off()
}

for (i in c("sB_Exclude3PCs","sI_Exclude3PCs","sB_Include3PCs","sI_Include3PCs")){
	dat <- get(i)
	n   <- nrow(dat)
  write.csv(dat, file=paste("./StateseQTL_",i, "_FDR.05_PEER10_", n, "_topSNPs_", d, ".csv", sep=""), row.names=FALSE)
}

## You have to re-read the *.csv files into R-studio otherwise the rD_Exclude3PCs$seqnames[i] will become factor and the following codes will fail
rD_Exclude3PCs <- read.csv("./ResponseQTL_rD_Exclude3PCs_FDR.05_PEER10_26_topSNPs_2016-09-23.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
rR_Exclude3PCs <- read.csv("./ResponseQTL_rR_Exclude3PCs_FDR.05_PEER10_22_topSNPs_2016-09-23.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
rD_Include3PCs <- read.csv("./ResponseQTL_rD_Include3PCs_FDR.05_PEER10_20_topSNPs_2016-09-23.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
rR_Include3PCs <- read.csv("./ResponseQTL_rR_Include3PCs_FDR.05_PEER10_12_topSNPs_2016-09-23.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
sB_Exclude3PCs <- read.csv("./StateseQTL_sB_Exclude3PCs_FDR.05_PEER10_5470_topSNPs_2016-09-27.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
sI_Exclude3PCs <- read.csv("./StateseQTL_sI_Exclude3PCs_FDR.05_PEER10_5184_topSNPs_2016-09-27.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
sB_Include3PCs <- read.csv("./StateseQTL_sB_Include3PCs_FDR.05_PEER10_5422_topSNPs_2016-09-27.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
sI_Include3PCs <- read.csv("./StateseQTL_sI_Include3PCs_FDR.05_PEER10_5290_topSNPs_2016-09-27.csv", stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)

rD26Ex3PCsSE <- get(load(paste(sourceDir, "/ResponseQTL/Without3PCs/imputed_Delta/SummarizedExpObject_withoutgQTLstatsRegOut_Delta_imputed_PEER-10_2016-09-21.rda", sep = "")))
for (i in 1:nrow(rD_Exclude3PCs)){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Delta_eGene_", chrnum, "_", rD_Exclude3PCs$probeid[i], "-3PCs.pdf", sep=""))
  CHR=paste("chr", rD_Exclude3PCs$seqnames[i], sep="")
  if ( rD_Exclude3PCs$seqnames[i] < 10 ) {TF=paste("tf0",rD_Exclude3PCs$seqnames[i],sep="")} else {TF=paste("tf",rD_Exclude3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- rD26Ex3PCsSE[seqnames(rD26Ex3PCsSE) == CHR, ]
  pV <- as.character(rD_Exclude3PCs$pValue[i])
  eqBox2(gene=rD_Exclude3PCs$probeid[i], GRanges(as.character(rD_Exclude3PCs$seqnames[i]), 
         IRanges(rD_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(rD_Exclude3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(rD_Exclude3PCs$probeid[i]))
  GT <- eqDesc2(gene=rD_Exclude3PCs$probeid[i], GRanges(as.character(rD_Exclude3PCs$seqnames[i]), 
         IRanges(rD_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}

rR22Ex3PCsSE <- get(load(paste(sourceDir, "/ResponseQTL/Without3PCs/imputed_Ratio/SummarizedExpObject_withoutgQTLstatsRegOut_Ratio_imputed_PEER-10_2016-09-20.rda", sep = "")))
for (i in 1:nrow(rR_Exclude3PCs)){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Ratio_eGene_", chrnum, "_", rR_Exclude3PCs$probeid[i], "-3PCs.pdf", sep=""))
  CHR=paste("chr", rR_Exclude3PCs$seqnames[i], sep="")
  if ( rR_Exclude3PCs$seqnames[i] < 10 ) {TF=paste("tf0",rR_Exclude3PCs$seqnames[i],sep="")} else {TF=paste("tf",rR_Exclude3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- rR22Ex3PCsSE[seqnames(rR22Ex3PCsSE) == CHR, ]
  pV <- as.character(rR_Exclude3PCs$pValue[i])
  eqBox2(gene=rR_Exclude3PCs$probeid[i], GRanges(as.character(rR_Exclude3PCs$seqnames[i]), 
         IRanges(rR_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(rR_Exclude3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(rR_Exclude3PCs$probeid[i]))
  GT <- eqDesc2(gene=rR_Exclude3PCs$probeid[i], GRanges(as.character(rR_Exclude3PCs$seqnames[i]), 
         IRanges(rR_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}

rD20In3PCsSE <- get(load(paste(sourceDir, "/ResponseQTL/With3PCs/imputed_Delta/SummarizedExpObject_withoutgQTLstatsRegOut_Delta_imputed_PEER-10_2016-09-20.rda", sep = "")))
for (i in 1:nrow(rD_Include3PCs)){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Delta_eGene_", chrnum, "_", rD_Include3PCs$probeid[i], "+3PCs.pdf", sep=""))
  CHR=paste("chr", rD_Include3PCs$seqnames[i], sep="")
  if ( rD_Include3PCs$seqnames[i] < 10 ) {TF=paste("tf0",rD_Include3PCs$seqnames[i],sep="")} else {TF=paste("tf",rD_Include3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- rD20In3PCsSE[seqnames(rD20In3PCsSE) == CHR, ]
  pV <- as.character(rD_Include3PCs$pValue[i])
  eqBox2(gene=rD_Include3PCs$probeid[i], GRanges(as.character(rD_Include3PCs$seqnames[i]), 
         IRanges(rD_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(rD_Include3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(rD_Include3PCs$probeid[i]))
  GT <- eqDesc2(gene=rD_Include3PCs$probeid[i], GRanges(as.character(rD_Include3PCs$seqnames[i]), 
         IRanges(rD_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}

rR12In3PCsSE <- get(load(paste(sourceDir, "/ResponseQTL/With3PCs/imputed_Ratio/SummarizedExpObject_withoutgQTLstatsRegOut_Ratio_imputed_PEER-10_2016-09-20.rda", sep = "")))
for (i in 1:nrow(rR_Include3PCs)){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Ratio_eGene_", chrnum, "_", rR_Include3PCs$probeid[i], "+3PCs.pdf", sep=""))
  CHR=paste("chr", rR_Include3PCs$seqnames[i], sep="")
  if ( rR_Include3PCs$seqnames[i] < 10 ) {TF=paste("tf0",rR_Include3PCs$seqnames[i],sep="")} else {TF=paste("tf",rR_Include3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- rR12In3PCsSE[seqnames(rR12In3PCsSE) == CHR, ]
  pV <- as.character(rR_Include3PCs$pValue[i])
  eqBox2(gene=rR_Include3PCs$probeid[i], GRanges(as.character(rR_Include3PCs$seqnames[i]), 
         IRanges(rR_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(rR_Include3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(rR_Include3PCs$probeid[i]))
  GT <- eqDesc2(gene=rR_Include3PCs$probeid[i], GRanges(as.character(rR_Include3PCs$seqnames[i]), 
         IRanges(rR_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}

sB5470Ex3PCsSE <- get(load(paste(sourceDir, "/StateseQTL/Without3PCs/imputed_Baseline/SummarizedExpObject_withoutgQTLstatsRegOut_Baseline_imputed_PEER-10_2016-09-19.rda", sep = "")))
for (i in 1:50){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Baseline_eGene_", chrnum, "_", sB_Exclude3PCs$probeid[i], "-3PCs.pdf", sep=""))
  CHR=paste("chr", sB_Exclude3PCs$seqnames[i], sep="")
  if ( sB_Exclude3PCs$seqnames[i] < 10 ) {TF=paste("tf0",sB_Exclude3PCs$seqnames[i],sep="")} else {TF=paste("tf",sB_Exclude3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- sB5470Ex3PCsSE[seqnames(sB5470Ex3PCsSE) == CHR, ]
  pV <- as.character(sB_Exclude3PCs$pValue[i])
  eqBox2(gene=sB_Exclude3PCs$probeid[i], GRanges(as.character(sB_Exclude3PCs$seqnames[i]), 
         IRanges(sB_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(sB_Exclude3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(sB_Exclude3PCs$probeid[i]))
  GT <- eqDesc2(gene=sB_Exclude3PCs$probeid[i], GRanges(as.character(sB_Exclude3PCs$seqnames[i]), 
         IRanges(sB_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}

sI5184Ex3PCsSE <- get(load(paste(sourceDir, "/StateseQTL/Without3PCs/imputed_Ischemia/SummarizedExpObject_withoutgQTLstatsRegOut_Ischemia_imputed_PEER-10_2016-09-19.rda", sep = "")))
for (i in 1:50){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Ischemia_eGene_", chrnum, "_", sI_Exclude3PCs$probeid[i], "-3PCs.pdf", sep=""))
  CHR=paste("chr", sI_Exclude3PCs$seqnames[i], sep="")
  if ( sI_Exclude3PCs$seqnames[i] < 10 ) {TF=paste("tf0",sI_Exclude3PCs$seqnames[i],sep="")} else {TF=paste("tf",sI_Exclude3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- sI5184Ex3PCsSE[seqnames(sI5184Ex3PCsSE) == CHR, ]
  pV <- as.character(sI_Exclude3PCs$pValue[i])
  eqBox2(gene=sI_Exclude3PCs$probeid[i], GRanges(as.character(sI_Exclude3PCs$seqnames[i]), 
         IRanges(sI_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(sI_Exclude3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(sI_Exclude3PCs$probeid[i]))
  GT <- eqDesc2(gene=sI_Exclude3PCs$probeid[i], GRanges(as.character(sI_Exclude3PCs$seqnames[i]), 
         IRanges(sI_Exclude3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}
                          
sB5422In3PCsSE <- get(load(paste(sourceDir, "/StateseQTL/With3PCs/imputed_Baseline/SummarizedExpObject_withoutgQTLstatsRegOut_Baseline_imputed_PEER-10_2016-09-19.rda", sep = "")))
for (i in 1:50){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Baseline_eGene_", chrnum, "_", sB_Include3PCs$probeid[i], "+3PCs.pdf", sep=""))
  CHR=paste("chr", sB_Include3PCs$seqnames[i], sep="")
  if ( sB_Include3PCs$seqnames[i] < 10 ) {TF=paste("tf0",sB_Include3PCs$seqnames[i],sep="")} else {TF=paste("tf",sB_Include3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- sB5422In3PCsSE[seqnames(sB5422In3PCsSE) == CHR, ]
  pV <- as.character(sB_Include3PCs$pValue[i])
  eqBox2(gene=sB_Include3PCs$probeid[i], GRanges(as.character(sB_Include3PCs$seqnames[i]), 
         IRanges(sB_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(sB_Include3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(sB_Include3PCs$probeid[i]))
  GT <- eqDesc2(gene=sB_Include3PCs$probeid[i], GRanges(as.character(sB_Include3PCs$seqnames[i]), 
         IRanges(sB_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}

sI5290In3PCsSE <- get(load(paste(sourceDir, "/StateseQTL/With3PCs/imputed_Ischemia/SummarizedExpObject_withoutgQTLstatsRegOut_Ischemia_imputed_PEER-10_2016-09-19.rda", sep = "")))
for (i in 1:50){
	if ( i < 10 ) {chrnum=paste("0",i,sep="")} else {chrnum=i}
  pdf(file=paste("Ischemia_eGene_", chrnum, "_", sI_Include3PCs$probeid[i], "+3PCs.pdf", sep=""))
  CHR=paste("chr", sI_Include3PCs$seqnames[i], sep="")
  if ( sI_Include3PCs$seqnames[i] < 10 ) {TF=paste("tf0",sI_Include3PCs$seqnames[i],sep="")} else {TF=paste("tf",sI_Include3PCs$seqnames[i],sep="")}
  TFR=get(TF)
  exprSE <- sI5290In3PCsSE[seqnames(sI5290In3PCsSE) == CHR, ]
  pV <- as.character(sI_Include3PCs$pValue[i])
  eqBox2(gene=sI_Include3PCs$probeid[i], GRanges(as.character(sI_Include3PCs$seqnames[i]), 
         IRanges(sI_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR,   xlab=paste(as.character(sI_Include3PCs$snp[i]), 
         "   pValue= ", pV, sep=""), ylab=as.character(sI_Include3PCs$probeid[i]))
  GT <- eqDesc2(gene=sI_Include3PCs$probeid[i], GRanges(as.character(sI_Include3PCs$seqnames[i]), 
         IRanges(sI_Include3PCs$start[i], width=1)), se=exprSE, tf=TFR)
  mtext(names(GT),side=3,line=2, at=c(1.75,2.00,2.25))
  mtext(c(GT[1],GT[2],GT[3]),side=3,line=1, at=c(1.75,2.00,2.25))
  dev.off()
}


timing <- proc.time()
message("# Finishing The entire R Script with running time (in unit of seconds) = ", timing, " seconds")
print(timing)
