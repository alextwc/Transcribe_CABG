#!/usr/bin/Rscript
message("############################################################################################################")
message("## R-script to run MatrixEQTL against the normalized (log(Exp+1.1)) State-specific RNA-Seq datasets       ##")
message("##                                                                                                        ##")
message("##  Author: Dr. Alex Tzuu-Wang Chang                                                                      ##")
message("##    Date: 2016-08-10                                                                                    ##")
message("##                                                                                                        ##")
message("## Version: 1.0 Using the Log2 transoformed RNA-Seq datasets instead of using the PEER-10 residuals       ##")
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
message("## Usage2 (for Channing Cluster): qbR -A "./LV135preX.csv ./LV133postX.csv" ./test.R                      ##")
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

require(MatrixEQTL)
d <- Sys.Date();
setwd("D:/BWHMS/Meetings02/2016-08-15");
base.dir <- "D:/BWHMS/Meetings02/2016-08-15/InputFiles_MatrixEQTL";
    odir <- getwd();
useModel <- modelLINEAR;

    ImputedSNP_file_name = paste(base.dir, "/Transcribe121_AllChr22_imputed_MAF0.15passed.recode.MEQTL.txt",  sep = "");
 expressionPre_file_name = paste(base.dir, "/exprPre.18214Genes.121Subjects.afterLog2.txt",  sep = "");
expressionPost_file_name = paste(base.dir, "/exprPost.18392Genes.121Subjects.afterLog2.txt", sep = "");
snps_location_file_name  = paste(base.dir, "/imputedSNPsPosition121.txt", sep = "");
gene_location_file_name  = paste(base.dir, "/genepos.txt", sep = "");
# output_file_name_cis   = tempfile();
# output_file_name_tra   = tempfile();
# covariates_file_name   = character(); # Set to character() for no covariates
 covariatesPre_file_name = paste(base.dir, "/121cvrtPre.txt",  sep = ""); 
covariatesPost_file_name = paste(base.dir, "/121cvrtPost.txt", sep = "");
output_file_name_cis     = paste(base.dir, "/MatrixEQTL.cis.imputed",   sep = "");
output_file_name_tra     = paste(base.dir, "/MatrixEQTL.trans.imputed", sep = "");
pvOutputThreshold_cis    = 1;
pvOutputThreshold_tra    = 0; # Don't run trans-eQTL analysis
# pvOutputThreshold_tra  = 1e-2;
snpspos                  = read.table(snps_location_file_name, header = FALSE, stringsAsFactors = FALSE);
genepos                  = read.table(gene_location_file_name, header = FALSE, stringsAsFactors = FALSE);
cisDist                  = 1e5;
errorCovariance          = numeric();

STATUS  <- c("Baseline", "Ischemia")
GeneExp <- c(expressionPre_file_name, expressionPost_file_name)
SNPs    <- c(ImputedSNP_file_name, ImputedSNP_file_name)
COVs    <- c(covariatesPre_file_name, covariatesPost_file_name)
rounds=seq(2)
for (i in rounds) {
  expression_file_name <- GeneExp[i]
  covariates_file_name <- COVs[i]
  SNP_file_name        <- SNPs[i]
  status               <- STATUS[i]
  
  message("###############################################################################")
  message("###############################################################################")
  message("################        Starting the round ", status, " now!       ################")
  message("###############################################################################")
  message("###############################################################################")
  print(Sys.time())
  message("")
  message("# Loading genotype data now #")
  snps                     = SlicedData$new();        # GxN matrix (G:SNPs , N:Subjects)
  snps$fileDelimiter       = "\t";                    # the Comma character
  snps$fileOmitCharacters  = "NA";                    # denote missing values
  snps$fileSkipRows        = 1;                       # one row of column labels
  snps$fileSkipColumns     = 1;                       # one column of row labels
  snps$fileSliceSize       = 50000;                   # read file in slices of 50,000 rows
  snps$LoadFile(SNP_file_name);                       # load SNPs for each 50,000 at a time
  message("")
  print(Sys.time())
  message("# filter out SNPs with MAF<=0.05 ...")
  # Note: here Minor allele is the minor one for the samples, not necessary the same one as the population)
  message("# filter out SNPs for those who has only one subject under minor allele homozygous")
  maf.list  = vector('list', length(snps))
  na.list   = vector('list', length(snps))
  mono.list = vector('list', length(snps))
  for(sl in 1:length(snps)) {
    slice = snps[[sl]];
    maf.list[[sl]]  = rowMeans(slice,na.rm=TRUE)/2; 
    maf.list[[sl]]  = pmin(maf.list[[sl]],1-maf.list[[sl]])
    mono.list[[sl]] = rowSums(slice==2, na.rm=TRUE)
    na.list[[sl]]   = is.na(rowMeans(slice));
  }
  maf  = unlist(maf.list)
  Mono = unlist(mono.list)
  na   = unlist(na.list)
  cat('SNPs before filtering:',nrow(snps), "\n")
  snps$RowReorder(Mono>1 & maf>0.05);
  cat('SNPs after filtering:',nrow(snps), "\n")
  rm(maf, Mono, na, maf.list, na.list, mono.list)
  print(Sys.time())
  message("")
  message("# Loading gene expression data now #")
  genes                    = SlicedData$new();        # GxN matrix (G:Geness , N:Subjects)
  genes$fileDelimiter      = "\t";                    # the Comma character
  genes$fileOmitCharacters = "NA";                    # denote missing values
  genes$fileSkipRows       = 1;                       # one row of column labels
  genes$fileSkipColumns    = 1;                       # one column of row labels
  genes$fileSliceSize      = 5000;                    # read file in slices of 5,000 rows
  genes$LoadFile(expression_file_name);               # load Genes expression for each 5,000 at a time
  message("")
  print(Sys.time())
  message("")
  message("# Loading covariates now #")
  cvrt                    = SlicedData$new();        # NxP matrix (N:Subjects , P:Phenotypes)
  cvrt$fileDelimiter      = "\t";                    # the Comma character
  cvrt$fileOmitCharacters = "NA";                    # denote missing values
  cvrt$fileSkipRows       = 1;                       # one row of column labels
  cvrt$fileSkipColumns    = 1;                       # one column of row labels
  cvrt$fileSliceSize      = 5000;                    # read file in slices of 5,000 rows
  if(length(covariates_file_name)>0) {               # load covariances now
    cvrt$LoadFile(covariates_file_name); }
  message("")
  print(Sys.time())
  message("")
  message("# Running MatrixEQTL on ", status, " dataset now #")
  print(Sys.time())
  ME1                     = Matrix_eQTL_main(
  snps                    = snps,
  gene                    = genes,
  cvrt                    = cvrt,
  output_file_name        = output_file_name_tra,
  pvOutputThreshold       = pvOutputThreshold_tra,
  useModel                = useModel,
  errorCovariance         = errorCovariance,
  verbose                 = FALSE,
  output_file_name.cis    = output_file_name_cis,
  pvOutputThreshold.cis   = pvOutputThreshold_cis,
  snpspos                 = snpspos,
  genepos                 = genepos,
  cisDist                 = cisDist,
  pvalue.hist             = TRUE);
  message("# Making the histogram of local and distant p-values #")
  pdf(file=paste("./pValueHistogram_imputed_", status, "_MatrixEQTL-RegressOut.pdf", sep=""), paper="usr", width = 0, height = 0)
   plot(ME1)
  dev.off()
  message(" ")
  print(Sys.time())
  ME2                      = Matrix_eQTL_main(
  snps                     = snps,
  gene                     = genes,
  cvrt                     = cvrt,
  output_file_name         = output_file_name_tra,
  pvOutputThreshold        = pvOutputThreshold_tra,
  useModel                 = useModel,
  errorCovariance          = errorCovariance,
  verbose                  = TRUE,
  output_file_name.cis     = paste0(output_file_name_cis, "_", status, "_MatrixEQTL-RegressOut.tab"),
  pvOutputThreshold.cis    = pvOutputThreshold_cis,
  snpspos                  = snpspos,
  genepos                  = genepos,
  cisDist                  = cisDist,
  pvalue.hist              = "qqplot",
  min.pv.by.genesnp        = FALSE,
  noFDRsaveMemory          = FALSE);
  # unlink(output_file_name_tra); # delete trans-eQTL result output_file
  # unlink(output_file_name_cis); # delete   cis-eQTL result output_file
  message(" ")
  message("# Making the Q-Q plot of local and distant p-values #")
  pdf(file=paste("./pValueQQplot_imputed_", status, "_MatrixEQTL-RegressOut.pdf", sep=""), paper="usr", width = 0, height = 0)
   plot(ME2)
  dev.off()
  message(" ")
  message("# Showing the results of MatrixEQTL analysis of ", status, " with covariates regressOut now #")
  cat('Analysis done in: ', ME2$time.in.sec, ' seconds', '\n');
  message("# The cis-eQTL analysis was done in ", (ME2$time.in.sec)/60, " minutes", "\n");
  message("# The number of total rows in cis-eQTL result is ", length(ME2$cis$eqtls$snps))
  message("# The number of harvested SNPs in cis-eQTL result is ", nlevels(factor(ME2$cis$eqtls$snps)))
  message("# The number of harvested eGenes in cis-eQTL result is ", nlevels(factor(ME2$cis$eqtls$gene)))
  message("# The cis-eQTL results are partially listing on the next lines: ")
  print(head(ME2$cis$eqtls))
  eQTL.FDR.05  <- nlevels(factor(ME2$cis$eqtls$snps[which(ME2$cis$eqtls$FDR<=0.05)]))
  eQTL.FDR.01  <- nlevels(factor(ME2$cis$eqtls$snps[which(ME2$cis$eqtls$FDR<=0.01)]))
  eGenes.FDR.05  <- nlevels(factor(ME2$cis$eqtls$gene[which(ME2$cis$eqtls$FDR<=0.05)]))
  eGenes.FDR.01  <- nlevels(factor(ME2$cis$eqtls$gene[which(ME2$cis$eqtls$FDR<=0.01)]))
  message("# The quantity of significant cis-eQTL with FDR smaller than 0.05 is ", eQTL.FDR.05)
  message("# The quantity of significant cis-eQTL with FDR smaller than 0.01 is ", eQTL.FDR.01)
  message("# The quantity of significant eGenes with FDR smaller than 0.05 is ",   eGenes.FDR.05)
  message("# The quantity of significant eGenes with FDR smaller than 0.01 is ",   eGenes.FDR.01)
  message(" ")
  message("# Running MatrixEQTL on ", status, " dataset without providing covariates now #")
  message(" ")
  print(Sys.time())
  message("# Reloading covariates with dummy file now #")
  dummyCvrt_file_name           = character();             # Set to character() for no covariates
  dummyCvrt                     = SlicedData$new();        # NxP matrix (N:Subjects , P:Phenotypes)
  dummyCvrt$fileDelimiter       = "\t";                    # the Comma character
  dummyCvrt$fileOmitCharacters  = "NA";                    # denote missing values
  dummyCvrt$fileSkipRows        = 1;                       # one row of column labels
  dummyCvrt$fileSkipColumns     = 1;                       # one column of row labels
  dummyCvrt$fileSliceSize       = 5000;                    # read file in slices of 5,000 rows
  if(length(dummyCvrt_file_name)>0) {                      # load covariances now
    dummyCvrt$LoadFile(dummyCvrt_file_name); }
  message("")
  print(Sys.time())
  message("")
  ME3                      = Matrix_eQTL_main(
  snps                     = snps,
  gene                     = genes,
  cvrt                     = dummyCvrt,
  output_file_name         = output_file_name_tra,
  pvOutputThreshold        = pvOutputThreshold_tra,
  useModel                 = useModel,
  errorCovariance          = errorCovariance,
  verbose                  = TRUE,
  output_file_name.cis     = paste0(output_file_name_cis, "_", status, "_MatrixEQTL-noRegressOut.tab"),
  pvOutputThreshold.cis    = pvOutputThreshold_cis,
  snpspos                  = snpspos,
  genepos                  = genepos,
  cisDist                  = cisDist,
  pvalue.hist              = "qqplot",
  min.pv.by.genesnp        = FALSE,
  noFDRsaveMemory          = FALSE);
  message(" ")
  message("# Making the Q-Q plot of local and distant p-values on NoRegressOut dataset now#")
   pdf(file=paste("./pValueQQplot_imputed_", status, "_MatrixEQTL-NoRegressOut.pdf", sep=""), paper="usr", width = 0, height = 0)
  plot(ME3)
  dev.off()
  message(" ")
  message("# Showing the results of MatrixEQTL analysis of ", status, " without covariates regressOut now #")
  cat('Analysis done in: ', ME3$time.in.sec, ' seconds', '\n');
  message("# The cis-eQTL analysis was done in ", (ME3$time.in.sec)/60, " minutes", "\n");
  message("# The number of total rows in cis-eQTL result is ", length(ME3$cis$eqtls$snps))
  message("# The number of harvested SNPs in cis-eQTL result is ", nlevels(factor(ME3$cis$eqtls$snps)))
  message("# The number of harvested eGenes in cis-eQTL result is ", nlevels(factor(ME3$cis$eqtls$gene)))
  message("# The cis-eQTL results are partially listing on the next lines: ")
  print(head(ME3$cis$eqtls))
  eQTL.FDR.05  <- nlevels(factor(ME3$cis$eqtls$snps[which(ME3$cis$eqtls$FDR<=0.05)]))
  eQTL.FDR.01  <- nlevels(factor(ME3$cis$eqtls$snps[which(ME3$cis$eqtls$FDR<=0.01)]))
  eGenes.FDR.05  <- nlevels(factor(ME3$cis$eqtls$gene[which(ME3$cis$eqtls$FDR<=0.05)]))
  eGenes.FDR.01  <- nlevels(factor(ME3$cis$eqtls$gene[which(ME3$cis$eqtls$FDR<=0.01)]))
  message("# The quantity of significant cis-eQTL with FDR smaller than 0.05 is ", eQTL.FDR.05)
  message("# The quantity of significant cis-eQTL with FDR smaller than 0.01 is ", eQTL.FDR.01)
  message("# The quantity of significant eGenes with FDR smaller than 0.05 is ",   eGenes.FDR.05)
  message("# The quantity of significant eGenes with FDR smaller than 0.01 is ",   eGenes.FDR.01)
  message(" ")
  print(Sys.time())
  message(" ")
}

print(Sys.time())

timing <- proc.time()
message("# Finishing The entire R Script with running time (in unit of seconds) = ", timing, " seconds")
print(timing)
