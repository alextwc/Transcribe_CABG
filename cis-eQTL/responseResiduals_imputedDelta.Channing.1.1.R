#!/usr/bin/Rscript
message("############################################################################################################")
message("## R-script to investigate and summarize the harvested Response eQTL eGenes from delta-residual files     ##")
message("## under different PEER-K factors (k=0,1,9,10,11,12,13,14,15) and finally to present in heatmap           ##")
message("##                                                                                                        ##")
message("##  Author: Dr. Alex Tzuu-Wang Chang                                                                      ##")
message("##    Date: 2016-01-26                                                                                    ##")
message("##                                                                                                        ##")
message("## Version: 1.1 Addiing the codes for ploting intra-conditional heatmap (contains 11 running conditions)  ##")
message("##             (noReg,ggReg,p00,p01,p09,p10,p11,p12,p13,p14,p15; also known as 2nd pass of PEER-K)        ##")
message("## Version: 1.0 delta-residuals from PEER-K(0,1,9,10,11,12,13,14,15)                                      ##")
message("##                                                                                                        ##")
message("## Require: R v3.2.1, python etc.                                                                         ##")
message("##                                                                                                        ##")
message("## When submit jobs at the outside of Channing:                                                           ##")
message("## (1) Login to capecod first (2) qsub -l lx yourScript.sh                                                ##")
message("## Switch to LINUX node: ssh -l retwc -o StrictHostKeyChecking=no nantucket.bwh.harvard.edu               ##")
message("##                                                                                                        ##")
message("## Run R in Channing:[ qs | (qrsh -l lx6,large) | (qrsh -l lx) ] then R; qrsh -l rstudio; R --vanilla;    ##")
message("## /udd/stvjc/VM/R-devel-dist/bin/R; ~stvjc/bin/Rdevel --vanilla;                                         ##")
message("## Run R script in Channing:                                                                              ##")
message("## nohup R --no-save < myRscript.R > myRscriptR.out &                                                     ##")
message("## R --vanilla < /udd/retwc/R/test.R > /udd/retwc/R/testR.out ; qsub -l lx ./runMyRscript.sh              ##")
message("##                                                                                                        ##")
message("## After qs or qrsh the R executable file is located at the following path:                               ##")
message("## /local/bin/R -> /app/R-3.2.0@i86-rhel6.0/bin/R                                                         ##")
message("## Rscript -> /app/R-3.2.0@i86-rhel6.0/bin/Rscript                                                        ##")
message("## Check Job status: qds retwc; qstat -r; qstat -f; qstat -f | grep retwc | less; qstat -ls; qhost;       ##")
message("##                                                                                                        ##")
message("## Usage1 (for Channing Cluster): ( Rscript ./test.R ./LV135preX.csv ./LV133postX.csv ) >& testR.log      ##")
message("## Usage2 (for Channing Cluster): qbR ./test.R ./LV135preX.csv ./LV133postX.csv                           ##")
message("## Result: Usage 1 works; Usage 2 (qbR command line) failed to pass the arguments                         ##")
message("############################################################################################################")
print(Sys.time())

R.version
R.home()
.libPaths()
sessionInfo()
require("pheatmap")

# setwd("D:/BWHMS/Meetings02/2016-01-18")
setwd("D:/BWHMS/Meetings02/2016-02-01")
PEERK <- c("00","01","09","10","11","12","13","14","15")
for(peerk in seq(9)) {
	knum <- PEERK[peerk]
	imputedDelta_noReg  <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/GGtoolsRegOut/ResponseQTL_deltaResidual_FDR.05_exprSE.csv",    sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_ggReg  <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/GGtoolsRegOut/ResponseQTL_deltaResidual_FDR.05_regExprSE.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p00Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-0.csv",  sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p01Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-1.csv",  sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p09Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-9.csv",  sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p10Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-10.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p11Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-11.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p12Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-12.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p13Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-13.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p14Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-14.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_p15Reg <- read.csv(file=paste("./residualResponse/cis-eQTL_results/PEER", knum, "_Residuals/deltaResidual_PEER-", knum, "/imputed/PEER_RegOut/ResponseQTL_deltaResidual_imputed_deltaResidual_PEER-", knum, "_FDR.05_PEER-15.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	imputedDelta_noReg  <-  imputedDelta_noReg[, c(16, 18, 11)]
	imputedDelta_ggReg  <-  imputedDelta_ggReg[, c(16, 18, 11)]
	imputedDelta_p00Reg <- imputedDelta_p00Reg[, c(16, 18, 11)]
	imputedDelta_p01Reg <- imputedDelta_p01Reg[, c(16, 18, 11)]
	imputedDelta_p09Reg <- imputedDelta_p09Reg[, c(16, 18, 11)]
	imputedDelta_p10Reg <- imputedDelta_p10Reg[, c(16, 18, 11)]
	imputedDelta_p11Reg <- imputedDelta_p11Reg[, c(16, 18, 11)]
	imputedDelta_p12Reg <- imputedDelta_p12Reg[, c(16, 18, 11)]
	imputedDelta_p13Reg <- imputedDelta_p13Reg[, c(16, 18, 11)]
	imputedDelta_p14Reg <- imputedDelta_p14Reg[, c(16, 18, 11)]
	imputedDelta_p15Reg <- imputedDelta_p15Reg[, c(16, 18, 11)]
	write.csv(imputedDelta_noReg,  file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_noReg.csv", sep=""),  row.names=FALSE)
	write.csv(imputedDelta_ggReg,  file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_ggReg.csv", sep=""),  row.names=FALSE)
	write.csv(imputedDelta_p00Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p00Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p01Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p01Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p09Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p09Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p10Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p10Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p11Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p11Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p12Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p12Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p13Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p13Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p14Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p14Reg.csv", sep=""), row.names=FALSE)
	write.csv(imputedDelta_p15Reg, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p15Reg.csv", sep=""), row.names=FALSE)
	
	test1 <- intersect(levels(factor(imputedDelta_noReg$probeid)),  levels(factor(imputedDelta_ggReg$probeid)))
	test2 <- intersect(levels(factor(imputedDelta_p00Reg$probeid)), levels(factor(imputedDelta_p01Reg$probeid)))
	test3 <- intersect(levels(factor(imputedDelta_p09Reg$probeid)), levels(factor(imputedDelta_p10Reg$probeid)))
	test4 <- intersect(levels(factor(imputedDelta_p11Reg$probeid)), levels(factor(imputedDelta_p12Reg$probeid)))
	test5 <- intersect(levels(factor(imputedDelta_p13Reg$probeid)), levels(factor(imputedDelta_p14Reg$probeid)))
	test1 <- intersect(levels(factor(imputedDelta_p15Reg$probeid)), test1)	
	test5 <- union(levels(factor(imputedDelta_noReg$probeid)),  levels(factor(imputedDelta_ggReg$probeid)))
	test6 <- union(levels(factor(imputedDelta_p00Reg$probeid)), levels(factor(imputedDelta_p01Reg$probeid)))
	test7 <- union(levels(factor(imputedDelta_p09Reg$probeid)), levels(factor(imputedDelta_p10Reg$probeid)))
	test8 <- union(levels(factor(imputedDelta_p11Reg$probeid)), levels(factor(imputedDelta_p12Reg$probeid)))
	test9 <- union(levels(factor(imputedDelta_p13Reg$probeid)), levels(factor(imputedDelta_p14Reg$probeid)))
	test5 <- union(levels(factor(imputedDelta_p15Reg$probeid)), test5)
	test2 <- intersect(test2, test3)
	test4 <- intersect(test4, test5)
	test2 <- intersect(test2, test4)
	test1 <- intersect(test2, test1)
	test6 <- union(test6, test7)
	test8 <- union(test8, test9)
	test9 <- union(test6, test8)
	test5 <- union(test9, test5)
	imputedDeltaResidual_eGene_intersectAll <- test1
	imputedDeltaResidual_eGene_unionAll     <- test5
	eGenes      <- c(imputedDeltaResidual_eGene_intersectAll, imputedDeltaResidual_eGene_unionAll)
	eGenes_From <- c(rep("IntersectionAll", length(imputedDeltaResidual_eGene_intersectAll)), rep("UnionAll", length(imputedDeltaResidual_eGene_unionAll)))
	eGeneDF     <- data.frame(eGenes_From=eGenes_From, eGenes=eGenes)
	write.csv(eGeneDF, file=paste("./imputedDelta_residualPEER", knum, "_eGenesTable-1.csv", sep=""), row.names=FALSE)
	iDI <- imputedDeltaResidual_eGene_intersectAll
	iDU <- imputedDeltaResidual_eGene_unionAll
	save(iDI, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_IntersectAll.rda", sep=""))
	save(iDU, file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_UnionAll.rda", sep=""))
	iDList  <- list(imputedDelta_noReg$probeid,imputedDelta_ggReg$probeid,imputedDelta_p00Reg$probeid,imputedDelta_p01Reg$probeid,imputedDelta_p09Reg$probeid,imputedDelta_p10Reg$probeid,imputedDelta_p11Reg$probeid,imputedDelta_p12Reg$probeid,imputedDelta_p13Reg$probeid,imputedDelta_p14Reg$probeid,imputedDelta_p15Reg$probeid)
	names(iDList) <- c("noReg","ggReg","Peer00","Peer01","Peer09","Peer10","Peer11","Peer12","Peer13","Peer14","Peer15")
	iDFrame <- data.frame(iDU,noReg=0,ggReg=0,p00Reg=0,p01Reg=0,p09Reg=0,p10Reg=0,p11Reg=0,p12Reg=0,p13Reg=0,p14Reg=0,p15Reg=0)
	rN <- nrow(iDFrame)
  for(k in seq(11)) {
     for(i in seq(rN)) {
	     if (iDFrame[i, 1] %in% iDList[[k]]) {
		     iDFrame[i, k+1] <- 1
	     } else {
		     iDFrame[i, k+1] <- 0
	     }
	   }
  }
  rownames(iDFrame) <- iDFrame[, 1]
  iDFrame <- iDFrame[, -1]
  iDFrame$Count <- rowSums(iDFrame==1)
  iDFrame <- iDFrame[order(iDFrame$Count, decreasing = TRUE), ]
  write.csv(iDFrame, file=paste("./imputedDelta_residualPEER", knum, "_eGenesTable-2.csv", sep=""), row.names=TRUE)
  iDFrame$Count <- NULL
  iDmatrix <- data.matrix(iDFrame)
  pdf(file=paste("./imputedDelta_residualPEER", knum, "_eGenesHeatMap.pdf", sep=""), paper="special", width = 9, height = 7)
   pheatmap(iDmatrix, legend=TRUE, legend_breaks=c(0,1), legend_labels=c("Insignificant","Significant"), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=FALSE, cellwidth=36, cellheight=5, main=paste("eGenes Heatmap of imputedDelta_residualPEER-", knum, " under 2nd Pass of PEER-K", sep=""), color=colorRampPalette(c("navy","white","firebrick3"))(4096), fontsize=10, fontsize_row=6, fontsize_col=10, fontsize_number=8)
  dev.off()  
}

iDU_PEER00 <- get(load("./InputFiles/imputedDelta_residualPEER00_UnionAll.rda"))
iDU_PEER01 <- get(load("./InputFiles/imputedDelta_residualPEER01_UnionAll.rda"))
iDU_PEER09 <- get(load("./InputFiles/imputedDelta_residualPEER09_UnionAll.rda"))
iDU_PEER10 <- get(load("./InputFiles/imputedDelta_residualPEER10_UnionAll.rda"))
iDU_PEER11 <- get(load("./InputFiles/imputedDelta_residualPEER11_UnionAll.rda"))
iDU_PEER12 <- get(load("./InputFiles/imputedDelta_residualPEER12_UnionAll.rda"))
iDU_PEER13 <- get(load("./InputFiles/imputedDelta_residualPEER13_UnionAll.rda"))
iDU_PEER14 <- get(load("./InputFiles/imputedDelta_residualPEER14_UnionAll.rda"))
iDU_PEER15 <- get(load("./InputFiles/imputedDelta_residualPEER15_UnionAll.rda"))
all_imputedDelta_Union <- list(iDU_PEER00,iDU_PEER01,iDU_PEER09,iDU_PEER10,iDU_PEER11,iDU_PEER12,iDU_PEER13,iDU_PEER14,iDU_PEER15)
names(all_imputedDelta_Union) <- c("iDUP00","iDUP01","iDUP09","iDUP10","iDUP11","iDUP12","iDUP13","iDUP14","iDUP15")
test1 <- union(all_imputedDelta_Union[[1]], all_imputedDelta_Union[[2]])
test2 <- union(all_imputedDelta_Union[[3]], all_imputedDelta_Union[[4]])
test3 <- union(all_imputedDelta_Union[[5]], all_imputedDelta_Union[[6]])
test4 <- union(all_imputedDelta_Union[[7]], all_imputedDelta_Union[[8]])
test1 <- union(test1, all_imputedDelta_Union[[9]])
test3 <- union(test3, test4)
test2 <- union(test1, test2)
test2 <- union(test2, test3)
iDeltaUnionAll <- test2
testdf <- data.frame(iDeltaUnionAll, iDPEER00_78=0,iDPEER01_38=0,iDPEER09_48=0,iDPEER10_45=0,iDPEER11_47=0,iDPEER12_48=0,iDPEER13_47=0,iDPEER14_45=0,iDPEER15_46=0)
rowN <- nrow(testdf)
for(k in seq(9)) {
   for(i in seq(rowN)) {
	   if (testdf[i, 1] %in% all_imputedDelta_Union[[k]]) {
		   testdf[i, k+1] <- 1
	   } else {
		   testdf[i, k+1] <- 0
	   }
	 }
}
rownames(testdf) <- testdf[, 1]
testdf <- testdf[, -1]
testdf$Count <- rowSums(testdf==1)
testdf <- testdf[order(testdf$Count, decreasing = TRUE), ]
write.csv(testdf, file=paste("./imputedDeltaResidualResult.csv", sep=""), row.names=TRUE)
testdf$Count <- NULL
testmatrix <- data.matrix(testdf)
pdf(file=paste("./imputedDeltaResidualsResonseGenes-1.pdf", sep=""), paper="special", width = 7, height = 9)
 pheatmap(testmatrix, legend=TRUE, legend_breaks=c(0,1), legend_labels=c("Insignificant","Significant"), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=FALSE, cellwidth=36, cellheight=5, main="Response eQTL genes by delta-residuals", color=colorRampPalette(c("navy","white","firebrick3"))(4096), fontsize=8, fontsize_row=6, fontsize_col=8, fontsize_number=6)
dev.off()
iD_noReg  <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_noReg.csv",  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_ggReg  <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_ggReg.csv",  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p00Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p00Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p01Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p01Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p09Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p09Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p10Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p10Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p11Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p11Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p12Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p12Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p13Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p13Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p14Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p14Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iD_p15Reg <- read.csv(file="./InputFiles/imputedDelta_residualPEER00_imputedDelta_p15Reg.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
iDListP00 <- list(iD_noReg$probeid,iD_ggReg$probeid,iD_p00Reg$probeid,iD_p01Reg$probeid,iD_p09Reg$probeid,iD_p10Reg$probeid,iD_p11Reg$probeid,iD_p12Reg$probeid,iD_p13Reg$probeid,iD_p14Reg$probeid,iD_p15Reg$probeid)
names(iDListP00) <- c("noReg_00","ggReg_00","p00Reg_00","p01Reg_00","p09Reg_00","p10Reg_00","p11Reg_00","p12Reg_00","p13Reg_00","p14Reg_00","p15Reg_00")
iDListAll <- iDListP00
PEERK <- c("01","09","10","11","12","13","14","15")
for(peerk in seq(8)) {
	knum <- PEERK[peerk]
	iD_noReg  <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_noReg.csv",  sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_ggReg  <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_ggReg.csv",  sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p00Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p00Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p01Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p01Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p09Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p09Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p10Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p10Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p11Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p11Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p12Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p12Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p13Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p13Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p14Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p14Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iD_p15Reg <- read.csv(file=paste("./InputFiles/imputedDelta_residualPEER", knum, "_imputedDelta_p15Reg.csv", sep=""), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	iDListK   <- list(iD_noReg$probeid,iD_ggReg$probeid,iD_p00Reg$probeid,iD_p01Reg$probeid,iD_p09Reg$probeid,iD_p10Reg$probeid,iD_p11Reg$probeid,iD_p12Reg$probeid,iD_p13Reg$probeid,iD_p14Reg$probeid,iD_p15Reg$probeid)
	names(iDListK) <- c(paste("noReg_",knum,sep=""),paste("ggReg_",knum,sep=""),paste("p00Reg_",knum,sep=""),paste("p01Reg_",knum,sep=""),paste("p09Reg_",knum,sep=""),paste("p10Reg_",knum,sep=""),paste("p11Reg_",knum,sep=""),paste("p12Reg_",knum,sep=""),paste("p13Reg_",knum,sep=""),paste("p14Reg_",knum,sep=""),paste("p15Reg_",knum,sep=""))
	iDListAll <- append(iDListAll,iDListK) #99 objects
}
# names(iDListAll) <- c("iDListP00","iDListP01","iDListP09","iDListP10","iDListP11","iDListP12","iDListP13","iDListP14","iDListP15")
testdf <- data.frame(iDeltaUnionAll,noReg_00=0,ggReg_00=0,p00Reg_00=0,p01Reg_00=0,p09Reg_00=0,p10Reg_00=0,p11Reg_00=0,p12Reg_00=0,p13Reg_00=0,p14Reg_00=0,p15Reg_00=0
                                   ,noReg_01=0,ggReg_01=0,p00Reg_01=0,p01Reg_01=0,p09Reg_01=0,p10Reg_01=0,p11Reg_01=0,p12Reg_01=0,p13Reg_01=0,p14Reg_01=0,p15Reg_01=0
                                   ,noReg_09=0,ggReg_09=0,p00Reg_09=0,p01Reg_09=0,p09Reg_09=0,p10Reg_09=0,p11Reg_09=0,p12Reg_09=0,p13Reg_09=0,p14Reg_09=0,p15Reg_09=0
                                   ,noReg_10=0,ggReg_10=0,p00Reg_10=0,p01Reg_10=0,p09Reg_10=0,p10Reg_10=0,p11Reg_10=0,p12Reg_10=0,p13Reg_10=0,p14Reg_10=0,p15Reg_10=0
                                   ,noReg_11=0,ggReg_11=0,p00Reg_11=0,p01Reg_11=0,p09Reg_11=0,p10Reg_11=0,p11Reg_11=0,p12Reg_11=0,p13Reg_11=0,p14Reg_11=0,p15Reg_11=0
                                   ,noReg_12=0,ggReg_12=0,p00Reg_12=0,p01Reg_12=0,p09Reg_12=0,p10Reg_12=0,p11Reg_12=0,p12Reg_12=0,p13Reg_12=0,p14Reg_12=0,p15Reg_12=0
                                   ,noReg_13=0,ggReg_13=0,p00Reg_13=0,p01Reg_13=0,p09Reg_13=0,p10Reg_13=0,p11Reg_13=0,p12Reg_13=0,p13Reg_13=0,p14Reg_13=0,p15Reg_13=0
                                   ,noReg_14=0,ggReg_14=0,p00Reg_14=0,p01Reg_14=0,p09Reg_14=0,p10Reg_14=0,p11Reg_14=0,p12Reg_14=0,p13Reg_14=0,p14Reg_14=0,p15Reg_14=0
                                   ,noReg_15=0,ggReg_15=0,p00Reg_15=0,p01Reg_15=0,p09Reg_15=0,p10Reg_15=0,p11Reg_15=0,p12Reg_15=0,p13Reg_15=0,p14Reg_15=0,p15Reg_15=0)
rowN <- nrow(testdf)
for(k in seq(99)) {
   for(i in seq(rowN)) {
	   if (testdf[i, 1] %in% iDListAll[[k]]) {
		   testdf[i, k+1] <- 1
	   } else {
		   testdf[i, k+1] <- 0
	   }
	 }
}
rownames(testdf) <- testdf[, 1]
testdf <- testdf[, -1]
testdf$Count <- rowSums(testdf==1)
testdf <- testdf[order(testdf$Count, decreasing = TRUE), ]
write.csv(testdf, file=paste("./imputedDeltaResidualResultAll99.csv", sep=""), row.names=TRUE)
testdf$Count <- NULL
testmatrix <- data.matrix(testdf)
pdf(file=paste("./imputedDeltaResidualsResonseGenesAll99.pdf", sep=""), paper="special", width = 25, height = 9)
 pheatmap(testmatrix, legend=TRUE, legend_breaks=c(0,1), legend_labels=c("Insignificant","Significant"), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=FALSE, cellwidth=16, cellheight=4, main="Response eQTL genes from All 99 normalized delta-residuals", color=colorRampPalette(c("navy","white","firebrick3"))(4096), fontsize=16, fontsize_row=5, fontsize_col=6, fontsize_number=5)
dev.off()

timing <- proc.time()
message("# Finishing The entire R Script with running time (in unit of seconds) listed in the following lines ")
print(timing)
print(timestamp())
