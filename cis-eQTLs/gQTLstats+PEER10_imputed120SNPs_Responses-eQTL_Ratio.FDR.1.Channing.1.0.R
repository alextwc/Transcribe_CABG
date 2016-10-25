#!/usr/bin/Rscript
message("############################################################################################################")
message("## R-script to run PEER10 & gQTLstats on Transcribe118 RNA-Seq for Response-Ratio cis-eQTL analyses       ##")
message("##                                                                                                        ##")
message("##  Author: Dr. Alex Tzuu-Wang Chang                                                                      ##")
message("##    Date: 2016-09-27                                                                                    ##")
message("##                                                                                                        ##")
message("## Version: 1.0 For re-running of 2016-09-08 Response cis-eQTL                                            ##")
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

message("###############################################################################")
message("# Environmental Variables and Working directory settings, R packages loading  #")
message("###############################################################################")
print(Sys.time())
setwd("/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160923reRUN01")
require("peer", lib.loc="/udd/retwc/R")
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
# require("tools")
listMarts()
base.dir <- getwd();
    odir <- getwd();
set.seed(12345)
imputed_all <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_AllChr22_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_01  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr01_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_02  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr02_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_03  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr03_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_04  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr04_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_05  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr05_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_06  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr06_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_07  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr07_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_08  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr08_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_09  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr09_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_10  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr10_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_11  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr11_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_12  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr12_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_13  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr13_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_14  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr14_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_15  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr15_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_16  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr16_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_17  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr17_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_18  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr18_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_19  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr19_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_20  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr20_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_21  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr21_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_22  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr22_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_23  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh23.vcf-IR.vcf.gz"
vcfpath_24  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh24.vcf-IR.vcf.gz"
vcfpath_25  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh25.vcf-IR.vcf.gz"
vcfpath_26  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh26.vcf-IR.vcf.gz"
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

message("###############################################################################")
message("#              Starting the command-line arguments inputs!                    #")
message("###############################################################################")
print(Sys.time())
args <- commandArgs(TRUE)
message("# The first command-line argument will be treated as the baseline(PRE) RNA-Seq. ")
message("# The second command-line argument will be treated as the ischemia(POST) RNA-Seq. ")
message("# The third command-line argument will be treated as the corresponded phenotypes dataset. ")
if (length(args)==0) {
  phenotypes <- read.csv("./InputFiles/CovariancesAll.csv", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  exprPre    <- read.csv("./InputFiles/LV135preX.csv",      header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  exprPost   <- read.csv("./InputFiles/LV133postX.csv",     header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  message("# No command-line arguments found, read the file directly from csv files #")
  # stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==3) {
  exprPre    <- read.csv(args[1], header=T, check.names=FALSE, stringsAsFactors=FALSE)
  exprPost   <- read.csv(args[2], header=T, check.names=FALSE, stringsAsFactors=FALSE)
  phenotypes <- read.csv(args[3], header=T, check.names=FALSE, stringsAsFactors=FALSE)
}
# system("ls x*")
# files <- system("ls x*",intern=T)
message("# Finished the command-line arguments inputs! ")
print(Sys.time())

message("                                                                               ")
message("                                                                               ")
message("###############################################################################")
message("# Preparing two states specific RNA-Seq expression datasets (still named as   #")
message("# PRE & POST) from both of the original PRE and POST RNA-Seq datasets.        #")
message("# Applying GTEx filtering & log2 transformation on the prepared datasets to   #")
message("# normalize the data distribution and to clean the improper data points such  #")
message("# as NA, NaN and Inf et cetera!                                               #")
message("###############################################################################")
message("# Preparing the RNA-Seq expression datasets .... removing postfixed X, removing gene coordinates information .....")
EXP <- list(df1=exprPre, df2=exprPost)
EXP <- lapply(EXP, function(df) {    
    df$ugene_id <- sub('X$', '', df$ugene_id) # remove the postfixed "X" protection
	  df <- df[,-c(2:5)] # remove columns of gene coordinates information
	  rownames(df) <- df[,1]
	  df <- df[, -1]
    df
})
exprPre  <- EXP$df1
exprPost <- EXP$df2
rm(EXP)
message("                                                                               ")
message("###############################################################################")
message("# Making the subset group of the RNA-Seq datset to the sample size you specified now   #")
message("###############################################################################")
# print(Sys.time())
# message("###############################################################################")
# message("# Restricting genes only from Ch20~Ch23                                       #")
# message("###############################################################################")
# exprPre  <- subset(exprPre, chr=="chr20" | chr=="chr21" | chr=="chr22" | chr=="chr23") 
# exprPost <- subset(exprPost, chr=="chr20" | chr=="chr21" | chr=="chr22" | chr=="chr23")
# exprPre  <- read.csv("./InputFiles/LV135preX.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
# exprPost <- read.csv("./InputFiles/LV133postX.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
# phenotypesFeMale <- subset(phenotypes, exclude==0 & Sex=="F", select=-c(AXC_time,ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
# rownames(phenotypesFeMale) <- phenotypesFeMale[, 1]
# phenotypesFeMale <- phenotypesFeMale[, -1]
# exprPre  <- exprPre[,  intersect(rownames(phenotypesFeMale), colnames(exprPre))]
# exprPost <- exprPost[, intersect(rownames(phenotypesFeMale), colnames(exprPost))]
# message("# Now the dimension of FEMALE PRE dataset is ", dim(exprPre)[1], " x ", dim(exprPre)[2])
# message("# Now the dimension of FEMALE POST dataset is ", dim(exprPost)[1], " x ", dim(exprPost)[2])
# phenotypes <- read.csv("./InputFiles/CovariancesAll.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)                                       
# GenoCov121 <- read.table("./InputFiles/Transcribe121_AllChr22_imputed_MAF0.15passed.ProPC.coord",header=TRUE,stringsAsFactors=FALSE)
mergedID   <- read.table("./InputFiles/geneid.mergedID",  header = TRUE, stringsAsFactors = FALSE)
GenoCov120 <- as.data.frame(fread("./InputFiles/Transcribe120_AllChr22_imputed_MAF0.15passed.ProPC.coord"))
rownames(GenoCov120) <- GenoCov120[, 1]
GenoCov120 <- GenoCov120[, -1]
rownames(GenoCov120)=gsub("(B.*V)", "\\1_\\1", rownames(GenoCov120))
GenoCov120 <- GenoCov120[, -c(1:4)]
GenoCov120 <- GenoCov120[,  c(1:3)]
exprPre  <-  exprPre[, intersect(rownames(GenoCov120), colnames(exprPre))]  # 23079x120
exprPost <- exprPost[, intersect(rownames(GenoCov120), colnames(exprPost))] # 23079x120

message("###############################################################################")
message("# Here are the object attributions of the original PRE expression dataset")
print(mode(exprPre)); print(class(exprPre)); print(dim(exprPre)); print(length(exprPre)); 
message("# Here are the object attributions of the original POST expression dataset")
print(mode(exprPost)); print(class(exprPost)); print(dim(exprPost)); print(length(exprPost)); 
message("###############################################################################")
print(Sys.time())
message("# Starting the filtration of zero expression genes now ")
exprPreZero  <-   exprPre[(rowSums(exprPre==0.0)  == ncol(exprPre)),]
exprPostZero <- exprPost[(rowSums(exprPost==0.0) == ncol(exprPost)),]
message("# The quantity of the genes with zero expression across entire ", dim(exprPre)[2],  " subjects of PRE dataset is ", dim(exprPreZero)[1])
message("# The quantity of the genes with zero expression across entire ", dim(exprPost)[2], " subjects of POST dataset is ", dim(exprPostZero)[1])
# exprPre  <-  exprPre[setdiff(rownames(exprPre),  rownames(exprPreZero)),] 
# exprPost <- exprPost[setdiff(rownames(exprPost), rownames(exprPostZero)),]
message("# Removing those genes with zero expression now! ")
exprPre  <-  exprPre[!(rowSums(exprPre==0.0)  == ncol(exprPre)),]  # 21962x120
exprPost <- exprPost[!(rowSums(exprPost==0.0) == ncol(exprPost)),] # 22047x120
message("# Removing those zero expression genes out from PRE dataset now ")
message("# Removing those zero expression genes out from POST dataset now ")
message("# The quantity of the genes left in PRE dataset after trimmed off zero expression genes is ", dim(exprPre)[1])
message("# The quantity of the genes left in POST dataset after trimmed off zero expression genes is ", dim(exprPost)[1])
message("# Removing those genes with individual expression larger than 0.1 RPKM but expressed in less than 10 subjects of PRE group now ") #GTEx Threshold
message("# Removing those genes with individual expression larger than 0.1 RPKM but expressed in less than 10 subjects of POST group now ") #GTEx Threshold
exprPre  <-  exprPre[rowSums(exprPre>0.1)  > 10,] # 18212 genes left
exprPost <- exprPost[rowSums(exprPost>0.1) > 10,] # 18377 genes left
message("# The quantity of the genes left in PRE dataset after GTEx [#subjects(0.1 RPKM)>10] threshold filtration is ", dim(exprPre)[1])
message("# The quantity of the genes left in POST dataset after GTEx [#subjects(0.1 RPKM)>10] threshold filtration is ", dim(exprPost)[1])
message("# Now it is ready to perform log2 transformation [log2(1.1+RPKM)] to normalize the datasets! ")
# message("# Now the dimension of PRE dataset is ", as.data.frame(print(dim(exprPre))))
# message("# Now the dimension of POST dataset is ", as.data.frame(print(dim(exprPost))))
message("# Now the dimension of PRE dataset is ", dim(exprPre)[1], " x ", dim(exprPre)[2])
message("# Now the dimension of POST dataset is ", dim(exprPost)[1], " x ", dim(exprPost)[2])
rm(exprPreZero)
rm(exprPostZero)
print(Sys.time())

# message("###############################################################################")
# message("# Are the genes of PRE dataset the same with POST dataset? ", setequal(rownames(exprPre), rownames(exprPost)))
# message("# Are the subjects of PRE dataset the same with POST dataset? ", setequal(colnames(exprPre), colnames(exprPost)))
message("###############################################################################")
message("##           Exploring the data point quality and distribution now           ##")
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) existed in PRE dataset is ", sum(is.na(exprPre)))
message("# The amount of NaN (Not a Number; 0/0) existed in PRE dataset is ", sum(is.nan(as.matrix(exprPre))))
message("# The amount of Inf (Infinity; caused by N/0) existed in PRE dataset is ", sum(is.infinite(as.matrix(exprPre))))
message("# The amount of zero in PRE dataset is ", sum(exprPre==0, na.rm=TRUE))
message("# The amount of negative values in PRE dataset is ", sum(exprPre < 0, na.rm=TRUE))
message("# The amount of positive values in PRE dataset is ", sum(exprPre > 0, na.rm=TRUE))
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) existed in POST dataset is ", sum(is.na(exprPost)))
message("# The amount of NaN (Not a Number; 0/0) existed in POST dataset is ", sum(is.nan(as.matrix(exprPost))))
message("# The amount of Inf (Infinity; caused by N/0) existed in POST dataset is ", sum(is.infinite(as.matrix(exprPost))))
message("# The amount of zero in POST dataset is ", sum(exprPost==0, na.rm=TRUE))
message("# The amount of negative values in POST dataset is ", sum(exprPost < 0, na.rm=TRUE))
message("# The amount of positive values in POST dataset is ", sum(exprPost > 0, na.rm=TRUE))
message("###############################################################################")
message("###############################################################################")
message("# Here are the object attributions of the PRE expression dataset after GTEx filtration ")
print(mode(exprPre)); print(class(exprPre)); print(dim(exprPre)); print(length(exprPre)); 
message("# Here are the object attributions of the POST expression dataset after GTEx filtration ")
print(mode(exprPost)); print(class(exprPost)); print(dim(exprPost)); print(length(exprPost));
message("###############################################################################")
message("###############################################################################")
message("# Now it is ready to perform log2 transformation [log2(1.1+RPKM)] to normalize the datasets! ")
print(Sys.time())
exprPre  <- log2(exprPre+1.1)
exprPost <- log2(exprPost+1.1)
message("# Finishing log2 Transformation to create Logarithmic expression datasets now! ")
print(Sys.time())
message("###############################################################################")
message("# Here are the object attributions of PRE dataset after log2 transformation ")
print(mode(exprPre)); print(class(exprPre)); print(dim(exprPre)); print(length(exprPre)); 
message("# Here are the object attributions of POST dataset after log2 transformation ")
print(mode(exprPost)); print(class(exprPost)); print(dim(exprPost)); print(length(exprPost)); 
message("###############################################################################")
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) in PRE dataset after log2 transformation is ", sum(is.na(exprPre)))
message("# The amount of NaN (Not a Number; 0/0) in PRE dataset after log2 transformation is ", sum(is.nan(as.matrix(exprPre))))
message("# The amount of Inf (Infinity; caused by N/0) in PRE dataset after log2 transformation is ", sum(is.infinite(as.matrix(exprPre))))
message("# The amount of zero in PRE dataset after log2 transformation is ", sum(exprPre==0, na.rm=TRUE))
message("# The amount of negative values in PRE dataset after log2 transformation is ", sum(exprPre < 0, na.rm=TRUE))
message("# The amount of positive values in PRE dataset after log2 transformation is ", sum(exprPre > 0, na.rm=TRUE))
message("# The Normalized PRE dataset has ", sum(exprPre==max(exprPre), na.rm=TRUE), " Max values and the Max value is ", max(exprPre, na.rm=TRUE))
message("# The Normalized PRE dataset has ", sum(exprPre==min(exprPre), na.rm=TRUE), " min values and the min value is ", min(exprPre, na.rm=TRUE))
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) in POST dataset after log2 transformation is ", sum(is.na(exprPost)))
message("# The amount of NaN (Not a Number; 0/0) in POST dataset after log2 transformation is ", sum(is.nan(as.matrix(exprPost))))
message("# The amount of Inf (Infinity; caused by N/0) in POST dataset after log2 transformation is ", sum(is.infinite(as.matrix(exprPost))))
message("# The amount of zero in POST dataset after log2 transformation is ", sum(exprPost==0, na.rm=TRUE))
message("# The amount of negative values in POST dataset after log2 transformation is ", sum(exprPost < 0, na.rm=TRUE))
message("# The amount of positive values in POST dataset after log2 transformation is ", sum(exprPost > 0, na.rm=TRUE))
message("# The Normalized POST dataset has ", sum(exprPost==max(exprPost), na.rm=TRUE), " Max values and the Max value is ", max(exprPost, na.rm=TRUE))
message("# The Normalized POST dataset has ", sum(exprPost==min(exprPost), na.rm=TRUE), " min values and the min value is ", min(exprPost, na.rm=TRUE))
message("###############################################################################")
message("###############################################################################")
message("# Saving the normalized datasets (by log2 transformation) after low-level processing now ")
# write.csv(NormalizedDelta, file=paste("./NormalizedDelta_", d, ".csv", sep=""), row.names=TRUE)
# d <- gsub("-", "_", Sys.Date())
print(Sys.time())
d <- Sys.Date()
geneQuanPre  <- dim(exprPre)[1]
subsQuanPre  <- dim(exprPre)[2]
geneQuanPost <- dim(exprPost)[1]
subsQuanPost <- dim(exprPost)[2]
write.csv(exprPre,  file=paste("./InputFiles/exprPre.",  geneQuanPre,  "Genes.",  subsQuanPre, "Subjects.afterLog2_", d, ".csv", sep=""), row.names=TRUE)
write.csv(exprPost, file=paste("./InputFiles/exprPost.", geneQuanPost, "Genes.", subsQuanPost, "Subjects.afterLog2_", d, ".csv", sep=""), row.names=TRUE)
rm(geneQuanPre)
rm(subsQuanPre)
rm(geneQuanPost)
rm(subsQuanPost)
print(Sys.time())
message("###############################################################################")
message("                                                                               ")
message("##################################################################################")
message("# Preparing Gene Coordinates and phenotypes for GRange & SummarizedExperiment    #")
message("# The regressOut(Sex, Age, seqCenter, reads, PC1~PC3) is prepared for gQTLstats  #")
message("##################################################################################")
print(Sys.time())
mergedID <- read.table("./InputFiles/geneid.mergedID", header = TRUE, stringsAsFactors = FALSE)
mergedID$ugene_id <- make.names(mergedID$ugene_id, unique=TRUE)
# phenotypes <- read.csv("./InputFiles/Covariances.csv", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
phenotypesPre  <- subset(phenotypes, exclude==0, select=-c(AXC_time,ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
phenotypesPost <- subset(phenotypes, exclude==0, select=-c(AXC_time,ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
rownames(phenotypesPre)  <-  phenotypesPre[, 1]
rownames(phenotypesPost) <- phenotypesPost[, 1]
phenotypesPre  <-  phenotypesPre[,-1]
phenotypesPost <- phenotypesPost[,-1]
# phenotypesPre  <- merge(GenoCov120, phenotypesPre,  by="row.names", all.x=TRUE)
# phenotypesPost <- merge(GenoCov120, phenotypesPost, by="row.names", all.x=TRUE)
# rownames(phenotypesPre)  <-  phenotypesPre[, 1]
# rownames(phenotypesPost) <- phenotypesPost[, 1]
# phenotypesPre  <-  phenotypesPre[,-1]
# phenotypesPost <- phenotypesPost[,-1]
exprPre  <-  exprPre[, intersect(rownames(phenotypesPre),   colnames(exprPre))]
exprPost <- exprPost[, intersect(rownames(phenotypesPost), colnames(exprPost))]
phenotypesPre  <-  phenotypesPre[intersect(rownames(phenotypesPre),   colnames(exprPre)), ]
phenotypesPost <- phenotypesPost[intersect(rownames(phenotypesPost), colnames(exprPost)), ]
# message("# Now the dimension of PRE dataset after intersection with phenotypes is ", as.data.frame(print(dim(exprPre))))
# message("# Now the dimension of POST dataset after intersection with phenotypes is ", as.data.frame(print(dim(exprPost))))
message("# Now the dimension of PRE dataset after intersection with phenotypes is ", dim(exprPre)[1], " x ", dim(exprPre)[2])
message("# Now the dimension of POST dataset after intersection with phenotypes is ", dim(exprPost)[1], " x ", dim(exprPost)[2])
message("# Saving the prepared phenotypic datasets now ")
# save(phenotypesPre,  file="./InputFiles/phenotypesPre.rda")
# save(phenotypesPost, file="./InputFiles/phenotypesPost.rda")
# write.csv(phenotypesPre,  file=paste("./InputFiles/phenotypesPre_", d, ".csv", sep=""), row.names=TRUE)
# write.csv(phenotypesPost, file=paste("./InputFiles/phenotypesPost_", d, ".csv", sep=""), row.names=TRUE)
print(Sys.time())

message("                                                                               ")
message("                                                                               ")
message("###############################################################################")
message("# Preparing the covs.rda object by using reshape2 package to convert the long #")
message("# data format to wide data format. It is because PEER package can only accept #")
message("# the wide data format of covariance table                                    #")
message("###############################################################################")
print(Sys.time())
require(reshape2)
CovsPre  <- read.csv("./InputFiles/CovariancesAll.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
CovsPost <- read.csv("./InputFiles/CovariancesAll.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
CovsPre  = subset(CovsPre,  exclude==0, select=-c(AXC_time,ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
CovsPost = subset(CovsPost, exclude==0, select=-c(AXC_time,ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
categorical_varaibles = c("Sex", "SeqCenter", "Reads")
message("# The identified categorical variables include SEX, Sequencing-Center and Reads-length ")
for(x in categorical_varaibles) {CovsPre  = cbind(CovsPre,  value=1);  CovsPre[,x]=paste0(x,  CovsPre[,x]); CovsPre  = dcast(CovsPre,  as.formula(paste0("... ~ ", x)), fill=0);}
for(x in categorical_varaibles) {CovsPost = cbind(CovsPost, value=1); CovsPost[,x]=paste0(x, CovsPost[,x]); CovsPost = dcast(CovsPost, as.formula(paste0("... ~ ", x)), fill=0);}
rownames(CovsPre)  =  CovsPre[, 1]
rownames(CovsPost) = CovsPost[, 1]
CovsPre  <-  CovsPre[,-1]
CovsPost <- CovsPost[,-1]
# CovsPre  <- merge(GenoCov120, CovsPre,  by="row.names", all.x=TRUE)
# CovsPost <- merge(GenoCov120, CovsPost, by="row.names", all.x=TRUE)
# rownames(CovsPre)  =  CovsPre[, 1]
# rownames(CovsPost) = CovsPost[, 1]
# CovsPre  <-  CovsPre[,-1]
# CovsPost <- CovsPost[,-1]
CovsPre  <-   CovsPre[intersect(rownames(CovsPre),   colnames(exprPre)), ]
CovsPost <-  CovsPost[intersect(rownames(CovsPost),  colnames(exprPost)), ]
exprPre  <-  exprPre[, intersect(rownames(CovsPre),  colnames(exprPre))]
exprPost <- exprPost[, intersect(rownames(CovsPost), colnames(exprPost))]
# Covs subjects were reduced from 169 to 132; 37 persons were removed by intersection
# rownames(Covs)=gsub("(B.*V)", "\\1_\\1", rownames(Covs)) # change BxxxxV back to BxxxxV_BxxxxV
message("# Finished the intersection of subjects between Covs tables and expression datasets ")
message("# Now the Covs tables and expression datasets have the same subjects ")
message("# Finished making the Covs tables for PEER factor analysis ")
message("# Saving the prepared Covs objects now ")
# save(CovsPre,  file="./InputFiles/covsPre.rda")
# save(CovsPost, file="./InputFiles/covsPost.rda")
# write.csv(CovsPre,  file="./InputFiles/CovsPre.csv",  row.names=TRUE)
# write.csv(CovsPost, file="./InputFiles/CovsPost.csv", row.names=TRUE)
# print(Sys.time())
# if(file.exists(".RData")) load(".RData") else{
#     message("# Didn't find the RData file so please load RData ...")}

message("                                                                               ")
message("                                                                               ")
message("###############################################################################")
message("# Obtaining the residual files from the first-pass of PEER factor analyses by #")
message("# providing the known covariates plus 10 hidden covariates                    #")
message("###############################################################################")
print(Sys.time())

model = PEER()
PEER_setCovariates(model, as.matrix(CovsPre)) # Include known covariates (NxC matrix) into the model 
PEER_setPhenoMean(model, as.matrix(t(exprPre))) # Include the expression dataset t(GxN matrix) into the model; Set expression data.
print(dim(PEER_getPhenoMean(model)))
PEER_setNk(model, 10) # 135/4=33.75 , so I choose "n_unobserved_factors=35"
getNk <- PEER_getNk(model)
message("# Starting the iterations with PEER_getNk = ", getNk)
PEER_setNmax_iterations(model, 10000) # As default, PEER iterates through updates of every variable 1000 times, you can change this limitation.
PEER_update(model) # Starting training the model, observing convergence
factors = PEER_getX(model)
weights = PEER_getW(model)
getCovs = PEER_getCovariates(model) # Please disable it while running PEER without including the known covariates
precision = PEER_getAlpha(model)
residuals = t(PEER_getResiduals(model))  # convert to GxN
rownames(factors) = rownames(CovsPre)
rownames(getCovs) = rownames(CovsPre) # Please disable it while running PEER without including the known covariates
rownames(weights) = rownames(exprPre)
rownames(residuals) = rownames(exprPre)
colnames(residuals) = colnames(exprPre)
pdf(file="./diagnostics_1stPEER10_pre.pdf", paper="usr", width = 0, height = 0)
 PEER_plotModel(model) # Please disable it while running PEER without including the known covariates
 mtext(paste0("PEER factor = ", k=10), side=3, line=1)
dev.off()
d <- Sys.Date()
k=10
status="Baseline"
write.csv(residuals, file=paste("./residuals_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
write.csv(factors,   file=paste("./factors_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
# Please disable write.csv(getCovs, ....) while running PEER without including the known covariates
write.csv(getCovs,   file=paste("./getCovs_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
write.csv(weights,   file=paste("./weights_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
write.csv(precision,   file=paste("./precision_1stPEER-", k, "_", status, ".csv", sep=""), row.names=FALSE)
pdf(file=paste("./precision_1stPEER-", k, "_", status, ".pdf", sep=""))
 plot(precision)
 mtext(paste0("PEER factor = ", k), side=3, line=-1)
dev.off()
rm(model)
rm(factors)
rm(getCovs)
rm(precision)
rm(residuals)
rm(weights)

model = PEER()
PEER_setCovariates(model, as.matrix(CovsPost)) # Include known covariates (NxC matrix) into the model 
PEER_setPhenoMean(model, as.matrix(t(exprPost))) # Include the expression dataset t(GxN matrix) into the model; Set expression data.
print(dim(PEER_getPhenoMean(model)))
PEER_setNk(model, 10) # 135/4=33.75 , so I choose "n_unobserved_factors=35"
getNk <- PEER_getNk(model)
message("# Starting the iterations with PEER_getNk = ", getNk)
PEER_setNmax_iterations(model, 10000) # As default, PEER iterates through updates of every variable 1000 times, you can change this limitation.
PEER_update(model) # Starting training the model, observing convergence
factors = PEER_getX(model)
weights = PEER_getW(model)
getCovs = PEER_getCovariates(model) # Please disable it while running PEER without including the known covariates
precision = PEER_getAlpha(model)
residuals = t(PEER_getResiduals(model))  # convert to GxN
rownames(factors) = rownames(CovsPost)
rownames(getCovs) = rownames(CovsPost) # Please disable it while running PEER without including the known covariates
rownames(weights) = rownames(exprPost)
rownames(residuals) = rownames(exprPost)
colnames(residuals) = colnames(exprPost)
pdf(file="./diagnostics_1stPEER10_post.pdf", paper="usr", width = 0, height = 0)
 PEER_plotModel(model) # Please disable it while running PEER without including the known covariates
 mtext(paste0("PEER factor = ", k=10), side=3, line=1)
dev.off()
d <- Sys.Date()
k=10
status="Ischemia"
write.csv(residuals, file=paste("./residuals_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
write.csv(factors,   file=paste("./factors_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
# Please disable write.csv(getCovs, ....) while running PEER without including the known covariates
write.csv(getCovs,   file=paste("./getCovs_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
write.csv(weights,   file=paste("./weights_1stPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
write.csv(precision,   file=paste("./precision_1stPEER-", k, "_", status, ".csv", sep=""), row.names=FALSE)
pdf(file=paste("./precision_1stPEER-", k, "_", status, ".pdf", sep=""))
 plot(precision)
 mtext(paste0("PEER factor = ", k), side=3, line=-1)
dev.off()
rm(model)
rm(factors)
rm(getCovs)
rm(precision)
rm(residuals)
rm(weights)
message("###############################################################################")
message("# Finished the first PEER10 pass for Sex,Age,ReadLength,SeqCenter without     #")
message("# genoPC1~PC3                                                                 #")
message("# The acquired residual files will be processed by Delta and Ratio algorithm  #")
message("###############################################################################")
print(Sys.time())
message("                                                                               ")
message("                                                                               ")
message("###############################################################################")
message("# Making the Delta & Ratio residual files for the scond pass of PEER10        #")
message("###############################################################################")
print(Sys.time())
# BaselineResiduals_1stPEER10 <- read.csv("D:/BWHMS/Meetings02/2016-06-13/ResponseQTL/FirstPEER10/Baseline/residuals_PEER-10_Baseline.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# IschemiaResiduals_1stPEER10 <- read.csv("D:/BWHMS/Meetings02/2016-06-13/ResponseQTL/FirstPEER10/Ischemia/residuals_PEER-10_Ischemia.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
BaselineResiduals_1stPEER10 <- read.csv("./residuals_1stPEER-10_Baseline.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
IschemiaResiduals_1stPEER10 <- read.csv("./residuals_1stPEER-10_Ischemia.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
baselineResidual <- BaselineResiduals_1stPEER10
ischemiaResidual <- IschemiaResiduals_1stPEER10
print(Sys.time())
rownames(baselineResidual) <- baselineResidual[, 1]
rownames(ischemiaResidual) <- ischemiaResidual[, 1]
baselineResidual <- baselineResidual[, -1]
ischemiaResidual <- ischemiaResidual[, -1]
baselineResidual <- baselineResidual[, intersect(colnames(baselineResidual), colnames(ischemiaResidual))] 
baselineResidual <- baselineResidual[intersect(rownames(baselineResidual), rownames(ischemiaResidual)), ] 
ischemiaResidual <- ischemiaResidual[intersect(rownames(baselineResidual), rownames(ischemiaResidual)), ] 
message("# Finished the intersection between Baseline-residual and Ischemia-residual ")
message("# Now BaselineResiduals_1stPEER10 has ", dim(baselineResidual)[2], " subjects and IschemiaResiduals_1stPEER10 has ", dim(ischemiaResidual)[2], " subjects as well ")
message("# Now BaselineResiduals_1stPEER10 has ", dim(baselineResidual)[1], " genes and IschemiaResiduals_1stPEER10 has ", dim(ischemiaResidual)[1], " genes as well ")
message("###############################################################################")
message("###############################################################################")
message("## Checking the dimension equality between Baseline-residual and Ischemia-residual")
message("###############################################################################")
message("# Are their genes the same? ",    setequal(rownames(baselineResidual), rownames(ischemiaResidual)))
message("# Are their subjects the same? ", setequal(colnames(baselineResidual), colnames(ischemiaResidual)))
message("###############################################################################")
message("## Checking the data point quality of Baseline-residual")
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) in baselineResidual is ", sum(is.na(baselineResidual)))
message("# The amount of NaN (Not a Number; 0/0) in baselineResidual is ", sum(is.nan(as.matrix(baselineResidual))))
message("# The amount of Inf (Infinity; caused by N/0) in baselineResidual is ", sum(is.infinite(as.matrix(baselineResidual))))
message("# The amount of zero in baselineResidual is ", sum(baselineResidual==0, na.rm=TRUE))
message("# The amount of negative values in baselineResidual is ", sum(baselineResidual < 0, na.rm=TRUE))
message("# The amount of positive values in baselineResidual is ", sum(baselineResidual > 0, na.rm=TRUE))
message("# This baselineResidual dataset has ", sum(baselineResidual==max(baselineResidual), na.rm=TRUE), " Max values and the Max value is ", max(baselineResidual, na.rm=TRUE))
message("# This baselineResidual dataset has ", sum(baselineResidual==min(baselineResidual), na.rm=TRUE), " min values and the min value is ", min(baselineResidual, na.rm=TRUE))
message("###############################################################################")
message("## Checking the data point quality of Ischemia-residual")
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) in ischemiaResidual is ", sum(is.na(ischemiaResidual)))
message("# The amount of NaN (Not a Number; 0/0) in ischemiaResidual is ", sum(is.nan(as.matrix(ischemiaResidual))))
message("# The amount of Inf (Infinity; caused by N/0) in ischemiaResidual is ", sum(is.infinite(as.matrix(ischemiaResidual))))
message("# The amount of zero in ischemiaResidual is ", sum(ischemiaResidual==0, na.rm=TRUE))
message("# The amount of negative values in ischemiaResidual is ", sum(ischemiaResidual < 0, na.rm=TRUE))
message("# The amount of positive values in ischemiaResidual is ", sum(ischemiaResidual > 0, na.rm=TRUE))
message("# This ischemiaResidual dataset has ", sum(ischemiaResidual==max(ischemiaResidual), na.rm=TRUE), " Max values and the Max value is ", max(ischemiaResidual, na.rm=TRUE))
message("# This ischemiaResidual dataset has ", sum(ischemiaResidual==min(ischemiaResidual), na.rm=TRUE), " min values and the min value is ", min(ischemiaResidual, na.rm=TRUE))
message("###############################################################################")    
message("###############################################################################")
message("")
message("# Start making the delta-residual now ")
print(Sys.time())
message("# Adding fudge factor 20 to proceed the coordinate translation to eliminate the negative data points")
baselineResidual <- baselineResidual+20
ischemiaResidual <- ischemiaResidual+20
deltaResidual    <- (ischemiaResidual-baselineResidual)/((ischemiaResidual+baselineResidual)/2)
message("# Finished making the delta-residual now ")
message("###############################################################################")
message("## Exploring the data point quality and distribution of delta-residual now   ##")
message("###############################################################################")
message("# The following checkings are for delta-residual ")
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) in this deltaResidual is ", sum(is.na(deltaResidual)))
message("# The amount of NaN (Not a Number; 0/0) in this deltaResidual is ", sum(is.nan(as.matrix(deltaResidual))))
message("# The amount of Inf (Infinity; caused by N/0) in this deltaResidual is ", sum(is.infinite(as.matrix(deltaResidual))))
message("# The amount of zero in this deltaResidual is ", sum(deltaResidual==0, na.rm=TRUE))
message("# The amount of negative values in this deltaResidual is ", sum(deltaResidual < 0, na.rm=TRUE))
message("# The amount of positive values in this deltaResidual is ", sum(deltaResidual > 0, na.rm=TRUE))
message("# This deltaResidual dataset has ", sum(deltaResidual==max(deltaResidual), na.rm=TRUE), " Max values and the Max value is ", max(deltaResidual, na.rm=TRUE))
message("# This deltaResidual dataset has ", sum(deltaResidual==min(deltaResidual), na.rm=TRUE), " min values and the min value is ", min(deltaResidual, na.rm=TRUE))
message("###############################################################################")

if (sum(is.infinite(as.matrix(deltaResidual))) > 0 ) {
    	message("###############################################################################")
    	message("# Since the existed \"Inf\" data points will cause the unconvergent iterations  #")
    	message("# in PEER-K analyses therefore we have to convert the Inf points to zero now  #")
    	deltaResidual[deltaResidual=="Inf"] <- 0
    	deltaResidual[deltaResidual=="-Inf"] <- 0 
    	message("# Finished Converting Inf to zero now                                         #")
    	message("# Re-checking the data points quality now                                     #")
    	message("###############################################################################")
    	message("# The amount of NA (Not Available; missing values) in this deltaResidual is ", sum(is.na(deltaResidual)))
    	message("# The amount of NaN (Not a Number; 0/0) in this deltaResidual is ", sum(is.nan(as.matrix(deltaResidual))))
    	message("# The amount of Inf (Infinity; caused by N/0) in this deltaResidual is ", sum(is.infinite(as.matrix(deltaResidual))))
    	message("# The amount of zero in this deltaResidual is ", sum(deltaResidual==0, na.rm=TRUE))
    	message("# The amount of negative values in this deltaResidual is ", sum(deltaResidual < 0, na.rm=TRUE))
    	message("# The amount of positive values in this deltaResidual is ", sum(deltaResidual > 0, na.rm=TRUE))
    	message("# This deltaResidual dataset has ", sum(deltaResidual==max(deltaResidual), na.rm=TRUE), " Max values and the Max value is ", max(deltaResidual, na.rm=TRUE))
    	message("# This deltaResidual dataset has ", sum(deltaResidual==min(deltaResidual), na.rm=TRUE), " min values and the min value is ", min(deltaResidual, na.rm=TRUE))
    	message("###############################################################################")
    } else {
    	message("###############################################################################")
    	message("# There is no detected \"Inf\" data points in this dataset                       ")
    	message("###############################################################################")
    }
message("###############################################################################")
message("# Saving the delta residual file for now!") # 18070x120
write.csv(deltaResidual, file=paste("./InputFiles/deltaResidual_", d, ".csv", sep=""), row.names=TRUE)
message("# Finished saving the delta residual file")
message("###############################################################################")
message("###############################################################################")
message("")
message("# Start making the ratio-residual now")
print(Sys.time())
ratioResidual <- (ischemiaResidual/baselineResidual)
ratioResidual <- log2(ratioResidual+1.1)
message("# Finished making the ratio-residual")
message("###############################################################################")
message("## Exploring the data point quality and distribution of ratio-residual now   ##")
message("###############################################################################")
message("# The following checkings are for ratio-residual")
message("###############################################################################")
message("# The amount of NA (Not Available; missing values) in this ratioResidual is ", sum(is.na(ratioResidual)))
message("# The amount of NaN (Not a Number; 0/0) in this ratioResidual is ", sum(is.nan(as.matrix(ratioResidual))))
message("# The amount of Inf (Infinity; caused by N/0) in this ratioResidual is ", sum(is.infinite(as.matrix(ratioResidual))))
message("# The amount of zero in this ratioResidual is ", sum(ratioResidual==0, na.rm=TRUE))
message("# The amount of negative values in this ratioResidual is ", sum(ratioResidual < 0, na.rm=TRUE))
message("# The amount of positive values in this ratioResidual is ", sum(ratioResidual > 0, na.rm=TRUE))
message("# This ratioResidual dataset has ", sum(ratioResidual==max(ratioResidual), na.rm=TRUE), " Max values and the Max value is ", max(ratioResidual, na.rm=TRUE))
message("# This ratioResidual dataset has ", sum(ratioResidual==min(ratioResidual), na.rm=TRUE), " min values and the min value is ", min(ratioResidual, na.rm=TRUE))
message("###############################################################################")
if (sum(is.infinite(as.matrix(ratioResidual))) > 0 ) {
    	message("###############################################################################")
    	message("# Since the existed \"Inf\" data points will cause the unconvergent iterations  #")
    	message("# in PEER-K analyses therefore we have to convert the Inf points to zero now  #")
    	ratioResidual[ratioResidual=="Inf"] <- 0
    	ratioResidual[ratioResidual=="-Inf"] <- 0 
    	message("# Finished Converting Inf to zero now                                         #")
    	message("# Re-checking the data points quality now                                     #")
    	message("###############################################################################")
    	message("# The amount of NA (Not Available; missing values) in this ratioResidual is ", sum(is.na(ratioResidual)))
    	message("# The amount of NaN (Not a Number; 0/0) in this ratioResidual is ", sum(is.nan(as.matrix(ratioResidual))))
    	message("# The amount of Inf (Infinity; caused by N/0) in this ratioResidual is ", sum(is.infinite(as.matrix(ratioResidual))))
    	message("# The amount of zero in this ratioResidual is ", sum(ratioResidual==0, na.rm=TRUE))
    	message("# The amount of negative values in this ratioResidual is ", sum(ratioResidual < 0, na.rm=TRUE))
    	message("# The amount of positive values in this ratioResidual is ", sum(ratioResidual > 0, na.rm=TRUE))
    	message("# This ratioResidual dataset has ", sum(ratioResidual==max(ratioResidual), na.rm=TRUE), " Max values and the Max value is ", max(ratioResidual, na.rm=TRUE))
    	message("# This ratioResidual dataset has ", sum(ratioResidual==min(ratioResidual), na.rm=TRUE), " min values and the min value is ", min(ratioResidual, na.rm=TRUE))
    	message("###############################################################################")
    } else {
    	message("###############################################################################")
    	message("# There is no detected \"Inf\" data points in this dataset                       ")
    	message("###############################################################################")
    }
message("###############################################################################")
message("# Saving the ratio residual file now!") # 18070x120
write.csv(ratioResidual, file=paste("./InputFiles/ratioResidual_", d, ".csv", sep=""), row.names=TRUE)
message("# Finished saving the ratio residual file")
message("###############################################################################")
message("###############################################################################")
message("")
# deltaResidual <- read.csv("./InputFiles/deltaResidual_PEER10.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# ratioResidual <- read.csv("./InputFiles/ratioResidual_PEER10.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# print(Sys.time())
# rownames(deltaResidual) <- deltaResidual[, 1]
# rownames(ratioResidual) <- ratioResidual[, 1]
# deltaResidual <- deltaResidual[, -1]
# ratioResidual <- ratioResidual[, -1]
# mergedID <- read.table("./InputFiles/geneid.mergedID", header = TRUE, stringsAsFactors = FALSE)
# mergedID$ugene_id <- make.names(mergedID$ugene_id, unique = TRUE)
phenotypesResidual <- subset(phenotypes, exclude==0, select=-c(ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
rownames(phenotypesResidual)  <-  phenotypesResidual[, 1]
phenotypesResidual <- phenotypesResidual[, -1]
deltaResidual <- deltaResidual[, intersect(rownames(phenotypesResidual), colnames(deltaResidual))]
ratioResidual <- ratioResidual[, intersect(rownames(phenotypesResidual), colnames(ratioResidual))]
phenotypesResidual  <-  phenotypesResidual[intersect(rownames(phenotypesResidual), colnames(ratioResidual)), ]
message("# Now the dimension of deltaResidual after intersection with phenotypes is ", dim(deltaResidual)[1], " x ", dim(deltaResidual)[2])
message("# Now the dimension of ratioResidual after intersection with phenotypes is ", dim(ratioResidual)[1], " x ", dim(ratioResidual)[2])
message("# Saving the prepared phenotypic datasets now!")
write.csv(phenotypesResidual, file = paste("./InputFiles/PhenotypesResiduals_PEER10_", d, ".csv", sep = ""), row.names = TRUE)
pheno <- phenotypesResidual
print(Sys.time())
Cvrt <- read.csv("./InputFiles/CovariancesAll.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
Cvrt <- subset(Cvrt, exclude == 0, select = -c(Sex,SeqCenter,Reads,Age,ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
message("# We only included \"AXC_time\" as the known covariate this time for PEER factor regress out ")
rownames(Cvrt) <- Cvrt[, 1]
# Cvrt <- Cvrt[, -1]
Cvrt <- Cvrt[intersect(rownames(Cvrt), colnames(ratioResidual)), ]
message("# Subjects of PEER_Known_Covariate_Table were reduced from 169 to 118; 51 persons were removed by intersection ")
write.csv(Cvrt, file = paste("./InputFiles/KnownCovariateFor_2ndPEER10_", d, ".csv", sep = ""), row.names = FALSE)

message("                                                                               ")
message("#################################################################################")
message("# Making the 1st loop for non-imputed and imputed SNPs to run gQTLstats         #")
message("# Making the 2nd loop to run gQTLstats sequentially on Delta and Ratio datasets #")
message("# Making the 3rd loop to run gQTLstats with PEER-K factor regress out residual  #")
message("#################################################################################")
print(Sys.time())
odir <- getwd()
message("# The root directory of running this script is at ", odir)
STATES  <- "ratioResidual"
STATUS  <- "Ratio"
 PHENO  <- "pheno"
   COV  <- "Cvrt"
  snps  <- "imputed"


      message("###############################################################################")
      message("#                      Using imputed SNPs for gQTLstats                       #")
      message("#                    Setting the paths for imputed SNPs now                   #")
      message("###############################################################################")
imputed_all <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_AllChr22_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_01  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr01_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_02  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr02_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_03  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr03_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_04  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr04_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_05  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr05_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_06  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr06_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_07  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr07_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_08  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr08_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_09  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr09_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_10  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr10_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_11  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr11_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_12  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr12_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_13  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr13_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_14  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr14_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_15  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr15_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_16  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr16_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_17  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr17_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_18  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr18_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_19  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr19_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_20  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr20_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_21  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr21_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_22  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20160818/VCF-IR/Transcribe120_Chr22_imputed_MAF0.15passed.vcf-IR.vcf.gz"
vcfpath_23  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh23.vcf-IR.vcf.gz"
vcfpath_24  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh24.vcf-IR.vcf.gz"
vcfpath_25  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh25.vcf-IR.vcf.gz"
vcfpath_26  <- "/proj/regeps/regep00/studies/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/non-imputed/VCF-IR/Merged142FinalReportCh26.vcf-IR.vcf.gz"
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
    
    system(sprintf("mkdir %s", snps))
    cudir=paste(odir, "/", snps, sep="")
    setwd(cudir)
    odir2 <- getwd()
    
    N=seq(2)
    
        fdname <- "ratioResidual"
        status <- "Ratio"
        createdFolder <- paste(snps, "_", status, sep="")
        system(sprintf("mkdir %s", createdFolder))
        currdir=paste(odir2, "/", createdFolder, sep="")
        setwd(currdir)
        # system("type > zz.log")
        expeer <- ratioResidual
        phenos <- phenotypesResidual
          Covs <- Cvrt
          
        message("###############################################################################")
        message("# The current working directory is at ", currdir)
        message("###############################################################################")
        message("                                                                               ")
        message("###############################################################################")
        message("# Starting running PEER factor analyese at the following path ", currdir)
        message("###############################################################################")
        message("# Run eQTL with PEER k factors plus the known covariance                      #")
        message("# 1. run PEER with covariates, e.g. PEER_setCovariates(model, as.matrix(Covs))#")
        message("# 2. get residuals from PEER, e.g. residuals = PEER_getResiduals(model),      #")
        message("# 3. use the acquired residuals as new expression input                       #")
        message("# for gQTLstats with or without regressOut(Age, Sex, SeqCenter, Reads-length) #")
        message("###############################################################################")
        message("# Run PEER with k factors now")
        # require(peer)
        # load("./covs.rda")
        print(Sys.time())
        k=10
        
            message("# Starting the PEER factors analyses with k = ", k)
            model = PEER()
            PEER_setCovariates(model, as.matrix(Covs)) # Include known covariates (NxC matrix) into the model 
            PEER_setPhenoMean(model, as.matrix(t(expeer))) # Include the expression dataset t(GxN matrix) into the model; Set expression data.
            print(dim(PEER_getPhenoMean(model)))
            # PEER_setPhenoVar(model, as.matrix(expeer_variance)) # Include the expression data uncertainty.
            # PEER_setSparsityPrior(model, as.matrix(prior)) # Set prior connectivity
            # If prior connectivity is specified, setting the number of unobserved factors is not needed.
            # The number of factors in the prior information matrix over-rides any previous specification
            PEER_setNk(model,k) # 135/4=33.75 , so I choose "n_unobserved_factors=35"
            # PEER_setNk(model, 0)  # no hidden factors
            # Note that unlike PCA-type models, the number of unobserved factors is not crucial when no
            # prior is specified because PEER uses automatic relevance determination to choose a
            # suitable effective number of factors. Hence, n_unobserved_factors needs only to be set to a
            # sufficiently large value (for technical details see Stegle et al.1). If no prior information on the
            # magnitude of confounding effects is available, we recommend using 25% of the number of
            # individuals contained in the study but no more than 100 factors.
            getNk <- PEER_getNk(model)
            message("# Starting the iterations with PEER_getNk = ", getNk)
            # PEER_setAdd_mean(model, TRUE)
            PEER_setNmax_iterations(model, 10000) # As default, PEER iterates through updates of every variable 1000 times, you can change this limitation.
            # PEER finishes if the increase in lower bound on the model evidence ceases to change, or the variance of the residuals has stabilised
            # Lower bound will increase along with the iteration cycles and residuals variance will decrease along with the iteration cycles.
            # PEER_setTolerance(model, 1) # 
            # PEER_setVarTolerance(model, 0.1) # 
            # In general you can keep the bound tolerance fairly high, but should keep the variation tolerance quite low compared to the variance of 
            # the expression matrix. If unsure, use the default values (bound=0.001, variance=0.00001).
            PEER_update(model) # Starting training the model, observing convergence
            # If the model is not converged after the default 1,000 iterations, and the variance of the
            # residuals keeps decreasing, choose a higher value of iterations, e.g., 10,000.
            # A total of 100 iterations should be sufficient to reach convergence on most data sets.
            # X = PEER_getX(model)  # now, the getX() factors should include known covariates + mean
            # W = PEER_getW(model)
            factors = PEER_getX(model)
            weights = PEER_getW(model)
            getCovs = PEER_getCovariates(model) # Please disable it while running PEER without including the known covariates
            # cor(PEER_getX(model)[,1], PEER_getCovariates(model)[,1])
            precision = PEER_getAlpha(model)
            residuals = t(PEER_getResiduals(model))  # convert to GxN
            # residuals = PEER_getResiduals(model) # the code that Patrick Evans used which did not transpose the matrix
            rownames(factors) = rownames(Covs)
            rownames(getCovs) = rownames(Covs) # Please disable it while running PEER without including the known covariates
            rownames(weights) = rownames(expeer)
            rownames(residuals) = rownames(expeer)
            colnames(residuals) = colnames(expeer)
            
            # The following PDF block has to be remarked when run PEER without known covariates
            # pdf(file=paste("./diagnostics_2ndPEER-", k, "_", status, ".pdf", sep=""), paper="usr", width = 0, height = 0)
            #  PEER_plotModel(model) # Please disable it while running PEER without including the known covariates
            #  mtext(paste0("PEER factor = ", k), side=3, line=1)
            # dev.off()
            
            # write.table(residuals_k, file="./residuals_PEER_k.LV135pre.xls", sep="\t")
            # write.table(factors_k, file="./factors_PEER_k.LV135pre.xls", sep="\t")
            # write.table(weights_k, file="./weights_PEER_k.LV135pre.xls", sep="\t")
            write.csv(residuals, file=paste("./residuals_2ndPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
            write.csv(factors,   file=paste("./factors_2ndPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
            # Please disable write.csv(getCovs, ....) while running PEER without including the known covariates
            write.csv(getCovs,   file=paste("./getCovs_2ndPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
            write.csv(weights,   file=paste("./weights_2ndPEER-", k, "_", status, ".csv", sep=""), row.names=TRUE)
            write.csv(precision,   file=paste("./precision_2ndPEER-", k, "_", status, ".csv", sep=""), row.names=FALSE)
            pdf(file=paste("./precision_2ndPEER-", k, "_", status, ".pdf", sep=""))
             plot(precision)
             mtext(paste0("PEER factor = ", k), side=3, line=-1)
            dev.off()
            # expr <- residuals + apply(exprChr20_23, 1, mean) 
            # add rowMean to each gene expression residual value; 1 means applying the function "mean" on each row; 2 means on each column
            expr <- residuals
            rm(model)
            rm(factors)
            rm(getCovs)
            rm(precision)
            rm(residuals)
            rm(weights)
            # OBJ <- ls()
            # message("# Here are the currently existed objects: "); print(OBJ)
            message("# Finishing the cycle K = ", k)
            print(Sys.time())
            message("                                                                               ")
            message("###############################################################################")
            message("###############################################################################")
            message("# Starting making SE(SummarizedExperiment) objects with or without regressOut()") 
            message("# We don't apply regressOut() during gQTLstats since the PEER has done already ")
            message("###############################################################################")
            message("###############################################################################")
            print(Sys.time())
            expr <- as.data.frame(expr)
            expr$ugene_id <- rownames(expr)
            names(expr)[names(expr)=="expr$ugene_id"] <- "ugene_id"
            expr <- merge(mergedID, expr, by="ugene_id")
            expr$chr <- gsub("chrX", "chr23", expr$chr)
            expr$chr <- gsub("chrY", "chr24", expr$chr)
            message("# Applying the GRanges function now ")
            message("# Integrating the gene coordinates and subject-phenotypic information ")
            message("# into SummarizedExperiment objects (SE) now                          ")
            message("# Adding the information of 'rownames & colnames' into prepared SE objects now ")
            rowD <- GRanges(seqnames=expr$chr, IRanges(expr$start, expr$end), strand=expr$strand)
            assay <- data.matrix(expr[, -c(1:5)])
            exprSE <- SummarizedExperiment(assay, rowRanges=rowD, colData=DataFrame(pheno))
            rownames(exprSE) = expr$ugene_id
            colnames(exprSE) = rownames(pheno)
            save(exprSE, file=paste("./SummarizedExpObject_withoutgQTLstatsRegOut_", status, "_", snps, "_PEER-", k, "_", d, ".rda", sep=""))
            # if( i != 2) {
            #   message("# Do not include the AXC_time as the known covariate for Ischemia state now ")
            #   # regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter)+AXC_time+factor(Reads))
            #   # regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter)+AXC_time)
            #   regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter))
            #   save(regExprSE, file=paste("./SummarizedExpObject_gQTLstatsRegOut_", status, "_", snps, "_PEER-", k, "_", d, ".rda", sep=""))
            # } else { 
    	      #   message("# Removing the AXC_time from the known covariate list for Baseline state now ")
            #   #  regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter)+factor(Reads))
            #   regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter))
            #   save(regExprSE, file=paste("./SummarizedExpObject_gQTLstatsRegOut_", status, "_", snps, "_PEER-", k, "_", d, ".rda", sep=""))
            # }
            # regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter)+AXC_time+Reads)
            # regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter)+AXC_time)
            # regExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter))
            # save(regExprSE, file=paste("./regExprSE_", k, ".rda", sep=""))
            message("# Finished making SE(SummarizedExperiment) objects and saved them as   ")
            message("# exprSE (without RegressOut) and regExprSE (with regressOut)          ")
            print(Sys.time())
            message("                                                                               ")
            message("###############################################################################")
            message("###############################################################################")
            message("# Run gQTLstats with Genotyping VCF files on exprSE objects now!              #")
            message("# expreSE is the dataset prepared without repetitively regressingOut those    #")
            message("# covariance (Age, Sex, SeqCenter, AXC_time) during making SE objects         #")
            message("###############################################################################")
            message("###############################################################################")
            # require(gQTLstats)
            print(Sys.time())
            set.seed(12345)
            message("# Setting the seed = 12345 now ")
            message("# Setting the MAF threshold as MAF >= 0.15 ")
            message("# Setting cis-radius as 100kb ")
            message("# Setting the permutation times as 3 to calculate plug-in FDR values ")
            message("# Setting the lower bound genotype frequency (lbgtf) as 0.023 to filter out the one-patient minor allele homozygous issue ")
            message("# Starting running gQTLstats package now ")
            
            message("###############################################################################")
            message("                                                                               ")
            message("# This is the cycle for the state of ", status, " and PEER-K= ", k, " and SNPs used is ", snps)
            message("###############################################################################")
            message("# The used expression dataset is ", fdname)
            message("# gQTLstats is running with exprSE(without regressOut) object now ")
            message("# Checking for universal heterozygous loci for exclusion")
            message("# Checking ......   ")
    ChrSE <- c("exprSE02", "exprSE03", "exprSE04", "exprSE05", "exprSE06", "exprSE07", "exprSE08", "exprSE09", "exprSE10", "exprSE11", "exprSE12", "exprSE13", "exprSE14", "exprSE15", "exprSE16", "exprSE17", "exprSE18", "exprSE19", "exprSE20", "exprSE21", "exprSE22", "exprSE23")
    ChrN  <- c("chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23")
    tf    <- c("tf02", "tf03", "tf04", "tf05", "tf06", "tf07", "tf08", "tf09", "tf10", "tf11", "tf12", "tf13", "tf14", "tf15", "tf16", "tf17", "tf18", "tf19", "tf20", "tf21", "tf22", "tf23")
    REGOUT <- c("exprSE", "regExprSE")

        message("# Starting splitting exprSE(SummarizedExperiment) object ", fdname, " to individual chromosomal SE objects ")
        message("                                                                               ")
        message("###############################################################################")
        exprSE01 = exprSE[ seqnames(exprSE) == "chr1", ]
        print(table(seqnames(exprSE01[1:10,])))
        seqlevelsStyle(exprSE01) = "NCBI"
        run01 = cisEsts(exprSE01,tf01,nperm=3,lbmaf=0.15,lbgtf=0.023,cisradius=100000)
        run01$piFDR = gQTLstats:::pifdr(run01$chisq, c(run01$permScore_1, run01$permScore_2, run01$permScore_3))
        message("The amount of tested SNP-GENE association pairs belong to chr1 is ", length(run01))
        message("The amount of successful-FDR-calculations on association pairs belong to chr1 is ", length(run01$piFDR))
        df_fdr.1 <- run01[which(run01$piFDR <= 1)]
        df_fdr.1 <- as(df_fdr.1, "data.frame")
        save(run01, file=paste("./run_chr1.exprSE_PEER-", k, "_", status, "_FDR.1.rda", sep=""))
        message("Finishing the gQTLstats-eQTL on Chromosome chr1 ")
        
        CHR=seq(22)
        for(chr in CHR) {
            ChromosomeNum <- ChrN[chr]
            EXPSE <- ChrSE[chr]
            TF <- get(tf[chr])
            EXPSE = exprSE[ seqnames(exprSE) == ChromosomeNum, ]
            print(table(seqnames(EXPSE[1:10,])))
            seqlevelsStyle(EXPSE) = "NCBI"
            runChr = cisEsts(EXPSE,TF,nperm=3,lbmaf=0.15,lbgtf=0.023,cisradius=100000)
            runChr$piFDR = gQTLstats:::pifdr(runChr$chisq, c(runChr$permScore_1, runChr$permScore_2, runChr$permScore_3))
            message("The amount of tested SNP-GENE association pairs belong to ", ChromosomeNum, " is ", length(runChr))
            message("The amount of successful-FDR-calculations on association pairs belong to ", ChromosomeNum, " is ", length(runChr$piFDR))
            save(runChr, file=paste("./run_", ChromosomeNum, ".exprSE_PEER-", k, "_", status, "_FDR.1.rda", sep=""))
            runChr_fdr.1 <- runChr[which(runChr$piFDR <= 1)]
            df.runChr_fdr.1 <- as(runChr_fdr.1, "data.frame")
            df_fdr.1 <- rbind(df_fdr.1, df.runChr_fdr.1)
            message("Finishing the gQTLstats-eQTL on Chromosome ", ChromosomeNum)
        } 
        message("                                                                               ")
        message("                                                                               ")
        message("###############################################################################")
        message("###############################################################################")
        message("# Summarizing the cis-eQTL result of this cycle with PEER factor = ", k)
        exprSeQTL_fdr.1  <- df_fdr.1
        exprSeQTL_fdr.1  <- exprSeQTL_fdr.1[order(exprSeQTL_fdr.1$piFDR), ]
        exprSeQTL_fdr.1  <- exprSeQTL_fdr.1[order(exprSeQTL_fdr.1$chisq, decreasing = TRUE), ]
        exprSeQTL_fdr.05 <- exprSeQTL_fdr.1[which(exprSeQTL_fdr.1$piFDR <= 0.05), ]
        exprSeQTL_fdr.05 <- exprSeQTL_fdr.05[order(exprSeQTL_fdr.05$piFDR), ]
        exprSeQTL_fdr.05 <- exprSeQTL_fdr.05[order(exprSeQTL_fdr.05$chisq, decreasing = TRUE), ]
        message("# The number of the columns in obtained cis-eQTL table of FDR.1 is ",  length(exprSeQTL_fdr.1))
        message("# The number of the columns in obtained cis-eQTL table of FDR.05 is ", length(exprSeQTL_fdr.05))
        message("# The dimension of the obtained cis-eQTL table of FDR.1 is ", dim(exprSeQTL_fdr.1)[1], " x ", dim(exprSeQTL_fdr.1)[2])
        message("# The dimension of the obtained cis-eQTL table of FDR.05 is ", dim(exprSeQTL_fdr.05)[1], " x ", dim(exprSeQTL_fdr.05)[2])

            assPairs.1  <- dim(exprSeQTL_fdr.1)
            assPairs.05 <- dim(exprSeQTL_fdr.05)
            eGenes.1  <- nlevels(factor(exprSeQTL_fdr.1$probeid))
            eGenes.05 <- nlevels(factor(exprSeQTL_fdr.05$probeid))
            if(eGenes.05 != 0) {
              exprSeQTL_fdr.1$eGenesCount  <- nlevels(factor(exprSeQTL_fdr.1$probeid))
              exprSeQTL_fdr.05$eGenesCount <- nlevels(factor(exprSeQTL_fdr.05$probeid))
            }
            save(exprSeQTL_fdr.1,  file=paste("./ResponseQTL_", status, "_", snps, "_", fdname, "_FDR.1_PEER-", k, ".rda", sep=""))
            save(exprSeQTL_fdr.05, file=paste("./ResponseQTL_", status, "_", snps, "_", fdname, "_FDR.05_PEER-", k, ".rda", sep=""))
            write.csv(exprSeQTL_fdr.1,  file=paste("./ResponseQTL_", status, "_", snps, "_", fdname, "_FDR.1_PEER-", k, ".csv", sep=""), row.names=FALSE)
            write.csv(exprSeQTL_fdr.05, file=paste("./ResponseQTL_", status, "_", snps, "_", fdname, "_FDR.05_PEER-", k, ".csv", sep=""), row.names=FALSE)
            # print(unlist(res))
            # test <- as.data.frame(lapply(res, function(x) { x }))
            # test <- t(sapply(res, function(x) { x }))
            # test <- data.frame(res)
            # dim(as.data.frame(res[[1]]))  # [1] 1985   20
            # dim(as.data.frame(res[[2]]))  # [1] 1289   20
            # dim(as.data.frame(res[[3]]))  # [1]  948   20
            # dim(exprSeQTL_fdr.1)          # [1] 4222   20
            # exprSeQTL_fdr.1 <- unlist(res)
            # exprSeQTL_fdr.1 <- as(exprSeQTL_fdr.1, "data.frame")
            message("###############################################################################")
            message("###############################################################################")
            message("# The current Response-specific dataset is ", status)
            message("# gQTLstats was run under the condition without repetitively regressingOut ")
            message("# The current PEER factor = ", k)
            message("# The number of SNP-Gene association pairs have FDR <= 100% is ",  assPairs.1[1])
            message("# The number of SNP-Gene association pairs have FDR <=  5% is ", assPairs.05[1])
            message("# The count of cis-eQTL genes (eGenes) of FDR <= 100% is ", eGenes.1)
            message("# The count of cis-eQTL genes (eGenes) of FDR <=   5% is ", eGenes.05)
            message("# Finishing gQTLstats on current expression dataset ", fdname, " with PEER-K= ", k)
            print(Sys.time())
            print(proc.time())
            message("###############################################################################")
            message("###############################################################################")
            message("                                                                               ")
            message("                                                                               ")

timing <- proc.time()
message("# Finishing The entire R Script with running time (in unit of seconds) = ", timing, " seconds")
print(timing)
