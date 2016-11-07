path("############################################################################################################")
path("## Usage1 (for Channing Cluster): ( Rscript ./test.R ./LV135preX.csv ./LV133postX.csv ) >& testR.log      ##")
path("## Usage2 (for Channing Cluster): qbR -A "./LV135preX.csv ./LV133postX.csv" ./test.R                      ##")
path("## Result: Usage 1 works; Usage 2 (qbR command line) is also working now after adding -A                  ##")
path("                                                                                                            ")
path("############################################################################################################")
print(Sys.time())

R.version
R.home()
.libPaths()
sessionInfo()
require(reshape2)
require(pheatmap)
d <- Sys.Date();
base.dir <- getwd();
load("D:/BWHMS/Meetings02/2016-04-19/InputFiles/eGeneList.rda")

GWAS_Catalog_Cardiovascular_File = paste(base.dir, "/InputFiles/GWAS_Catalog_Cardiovascular.csv", sep = "")
GWAS_Catalog_Cardiovascular <- read.csv(GWAS_Catalog_Cardiovascular_File,  stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)

imputedBaseline_eGenes_InGWASCatalogCardiovascular <- imputedBaseline_eGenes[(imputedBaseline_eGenes %in% GWAS_Catalog_Cardiovascular$REPORTED_GENE)]
imputedIschemia_eGenes_InGWASCatalogCardiovascular <- imputedIschemia_eGenes[(imputedIschemia_eGenes %in% GWAS_Catalog_Cardiovascular$REPORTED_GENE)]
TranscribeSignificantDE_InGWASCatalogCardiovascular <- TranscribeCuffdiffSignificant[(TranscribeCuffdiffSignificant$gene_id %in% GWAS_Catalog_Cardiovascular$REPORTED_GENE), ]
write.table(imputedBaseline_eGenes_InGWASCatalogCardiovascular,    file="./imputedBaseline_eGenes_InGWASCatalogCardiovascular.txt",    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(imputedIschemia_eGenes_InGWASCatalogCardiovascular,    file="./imputedIschemia_eGenes_InGWASCatalogCardiovascular.txt",    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(GWAS_Catalog_Cardiovascular$REPORTED_GENE, file="./GWAS_Catalog_Cardiovascular.txt",    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.csv(TranscribeSignificantDE_InGWASCatalogCardiovascular, file=paste("./TranscribeSignificantDE_InGWASCatalogCardiovascular_", d, ".csv", sep=""), row.names=TRUE)
write.table(rownames(TranscribeGSEA), file="./ORA_Total_19290Genes_List.txt",    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")




> (intersect(rownames(TranscribeCuffdiff_iBaselineDE), rownames(TranscribeSignificantDE_InGWASCatalogCardiovascular)))
 [1] "FAM89A"  "TSPAN15" "PI16"    "ICAM1"   "LDB2"    "SCN2B"   "FASTKD2" "AZIN1"  
 [9] "PLCL2"   "MUT"     "APOB"    "ADCY9"   "CHRM2"   "PHACTR1" "ODZ2"    "CFH"    
[17] "CELSR2" 
> (intersect(rownames(TranscribeCuffdiff_iIschemiaDE), rownames(TranscribeSignificantDE_InGWASCatalogCardiovascular)))
 [1] "NOV"     "FAM89A"  "TSPAN15" "PI16"    "LDB2"    "SCN2B"   "MICAL2"  "AZIN1"  
 [9] "APOB"    "CCDC141" "CHRM2"   "PHACTR1" "ODZ2"    "CD163"  
> 




LV135pre_File  = paste(base.dir, "/InputFiles/LV135preX.csv",  sep = "");
LV133post_File = paste(base.dir, "/InputFiles/LV133postX.csv", sep = "");
TranscribeCuffDiff_File = paste(base.dir, "/InputFiles/Transcribe.gene_exp.diff", sep = "")
TranscribeCuffDiffSignificant_File = paste(base.dir, "/InputFiles/Transcribe.gene_exp.diff.Significant", sep = "")
eGeneList.rda_File = paste(base.dir, "/eGeneList.rda", sep = "");
LV_PEER10_Fudge20Log2Transformed_File = paste(base.dir, "/InputFiles/TranscribeGSEA_PEER10_Fudge20Log2Transformed2016-04-04.csv", sep = "");
# MouseHeart_Exp <- read.table(Mouse_Heart_Exp_file_name, header = TRUE, stringsAsFactors = FALSE, sep = ",")
# HumanLV_Exp      <- read.table(Human_LV_Exp_file_name,    stringsAsFactors = FALSE, header = FALSE, sep = "\t")
load(eGeneList.rda_File)
exprPre  <- read.csv(LV135pre_File,  stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
exprPost <- read.csv(LV133post_File, stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
TranscribeCuffdiff <- read.table(TranscribeCuffDiff_File, stringsAsFactors = FALSE, header = TRUE)
# TranscribeCuffdiffSignificant <- read.table(TranscribeCuffDiffSignificant_File, stringsAsFactors = FALSE, header = TRUE)
TranscribeCuffdiffSignificant <- subset(TranscribeCuffdiff, significant=="yes")
# TranscribeCuffdiffSignificant <- TranscribeCuffdiffSignificant[order(TranscribeCuffdiffSignificant$log2.fold_change.), ]
TranscribeCuffdiffSignificant <- TranscribeCuffdiffSignificant[order(TranscribeCuffdiffSignificant$log2.fold_change. , decreasing = TRUE), ]
write.csv(TranscribeCuffdiffSignificant, file=paste("./TranscribeCuffdiffSignificant_", d, ".csv", sep=""), row.names=TRUE)
length(intersect(TranscribeCuffdiffSignificant$gene_id, imputedBaseline_eGenes)) # 391
LV_PEER10_Fudge20_Log2 <- read.csv(LV_PEER10_Fudge20Log2Transformed_File, stringsAsFactors = FALSE, header = TRUE,  check.names = FALSE)
rownames(LV_PEER10_Fudge20_Log2) <- LV_PEER10_Fudge20_Log2[, 1]
rownames(TranscribeCuffdiff) <- TranscribeCuffdiff[, 2]
LV_PEER10_Fudge20_Log2 <- LV_PEER10_Fudge20_Log2[, -1]
TranscribeCuffdiff <- TranscribeCuffdiff[, -1]
path("# Preparing the RNA-Seq expression datasets .... removing postfixed X, removing gene coordinates information .....")
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

# path("                                                                               ")
# path("###############################################################################")
# path("# Preparing the FEMALE subset group of states specific RNA-Seq datasets now   #")
# path("###############################################################################")
# phenotypesFeMale <- subset(phenotypes, exclude==0 & Sex=="F", select=-c(AXC_time,ID,Sequence,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
# rownames(phenotypesFeMale) <- phenotypesFeMale[, 1]
# phenotypesFeMale <- phenotypesFeMale[, -1]
# exprPre  <- exprPre[,  intersect(rownames(phenotypesFeMale), colnames(exprPre))]
# exprPost <- exprPost[, intersect(rownames(phenotypesFeMale), colnames(exprPost))]
# path("# Now the dimension of FEMALE PRE dataset is ", dim(exprPre)[1], " x ", dim(exprPre)[2])
# path("# Now the dimension of FEMALE POST dataset is ", dim(exprPost)[1], " x ", dim(exprPost)[2])
# path("                                                                               ")
# path("###############################################################################")
# path("# Changing subject names to attach with Pre or Post status to run GSEA        #")
# path("###############################################################################")
 colnames(exprPre)=gsub("(B.*V)_(B.*V)", "\\1_Pre",  colnames(exprPre)) # change BxxxxV_BxxxxV to BxxxxV_Pre , 135 subjects
colnames(exprPost)=gsub("(B.*V)_(B.*V)", "\\1_Post", colnames(exprPost)) # change BxxxxV_BxxxxV to BxxxxV_Post, 133 subjects
TranscribeGSEA <- merge(exprPre, exprPost, by="row.names") # 23079 genes x 269 columns 
rownames(TranscribeGSEA) <- TranscribeGSEA[, 1]
TranscribeGSEA <- TranscribeGSEA[, -1] # 23079 genes x 268 subjects # 18258 x 267 (for PEER-10 residuals)
path("###############################################################################")
path("# Here are the object attributions of the original TranscribeGSEA expression dataset")
print(mode(TranscribeGSEA)); print(class(TranscribeGSEA)); print(dim(TranscribeGSEA)); print(length(TranscribeGSEA)); 
path("###############################################################################")
print(Sys.time())
path("# Starting the filtration of zero expression genes now ")
TranscribeGSEAZero <- TranscribeGSEA[(rowSums(TranscribeGSEA==0.0)  == ncol(TranscribeGSEA)),]
path("# The quantity of the genes with zero expression across entire ", dim(TranscribeGSEA)[2],  " subjects of TranscribeGSEA dataset is ", dim(TranscribeGSEAZero)[1])
path("# Removing those genes with zero expression now! ")
TranscribeGSEA <- TranscribeGSEA[!(rowSums(TranscribeGSEA==0.0)  == ncol(TranscribeGSEA)),] 
path("# Removing those zero expression genes out from TranscribeGSEA dataset now ")
path("# The quantity of the genes left in TranscribeGSEA dataset after trimmed off zero expression genes is ", dim(TranscribeGSEA)[1]) # 22251 genes left
path("# Removing those genes with individual expression larger than 0.1 RPKM but expressed in less than 10 subjects of TranscribeGSEA group now ") #GTEx Threshold
TranscribeGSEA <- TranscribeGSEA[rowSums(TranscribeGSEA>0.1)  > 10,] # 19290 genes left # 18217 genes left (for PEER-10 residuals)
path("# The quantity of the genes left in TranscribeGSEA dataset after GTEx [#subjects(0.1 RPKM)>10] threshold filtration is ", dim(TranscribeGSEA)[1])
path("# Now the dimension of TranscribeGSEA dataset is ", dim(TranscribeGSEA)[1], " x ", dim(TranscribeGSEA)[2]) # 19290 x 268 # 18217 x 267 (for PEER-10 residuals)
rm(TranscribeGSEAZero)
print(Sys.time())

path("###############################################################################")
path("##           Exploring the data point quality and distribution now           ##")
path("###############################################################################")
path("# The amount of NA (Not Available; missing values) existed in TranscribeGSEA dataset is ", sum(is.na(TranscribeGSEA)))
path("# The amount of NaN (Not a Number; 0/0) existed in TranscribeGSEA dataset is ", sum(is.nan(as.matrix(TranscribeGSEA))))
path("# The amount of Inf (Infinity; caused by N/0) existed in TranscribeGSEA dataset is ", sum(is.infinite(as.matrix(TranscribeGSEA))))
path("# The amount of zero in TranscribeGSEA dataset is ", sum(TranscribeGSEA==0, na.rm=TRUE))
path("# The amount of negative values in TranscribeGSEA dataset is ", sum(TranscribeGSEA < 0, na.rm=TRUE))
path("# The amount of positive values in TranscribeGSEA dataset is ", sum(TranscribeGSEA > 0, na.rm=TRUE))
path("# The TranscribeGSEA dataset has ", sum(TranscribeGSEA==max(TranscribeGSEA), na.rm=TRUE), " Max values and the Max value is ", max(TranscribeGSEA, na.rm=TRUE))
path("# The TranscribeGSEA dataset has ", sum(TranscribeGSEA==min(TranscribeGSEA), na.rm=TRUE), " min values and the min value is ", min(TranscribeGSEA, na.rm=TRUE))
path("###############################################################################")
path("###############################################################################")
path("# Here are the object attributions of the TranscribeGSEA expression dataset after GTEx filtration ")
print(mode(TranscribeGSEA)); print(class(TranscribeGSEA)); print(dim(TranscribeGSEA)); print(length(TranscribeGSEA)); 
path("###############################################################################")
path("# Now it is ready to perform log2 transformation [log2(1.1+RPKM)] to normalize the datasets! ")
print(Sys.time())
TranscribeGSEA  <- log2(TranscribeGSEA+1.1)
path("# Finishing log2 Transformation to create Logarithmic expression datasets now! ")
print(Sys.time())
path("###############################################################################")
path("# Here are the object attributions of TranscribeGSEA dataset after log2 transformation ")
print(mode(TranscribeGSEA)); print(class(TranscribeGSEA)); print(dim(TranscribeGSEA)); print(length(TranscribeGSEA)); 
path("###############################################################################")
path("# The amount of NA (Not Available; missing values) in TranscribeGSEA dataset after log2 transformation is ", sum(is.na(TranscribeGSEA)))
path("# The amount of NaN (Not a Number; 0/0) in TranscribeGSEA dataset after log2 transformation is ", sum(is.nan(as.matrix(TranscribeGSEA))))
path("# The amount of Inf (Infinity; caused by N/0) in TranscribeGSEA dataset after log2 transformation is ", sum(is.infinite(as.matrix(TranscribeGSEA))))
path("# The amount of zero in TranscribeGSEA dataset after log2 transformation is ", sum(TranscribeGSEA==0, na.rm=TRUE))
path("# The amount of negative values in TranscribeGSEA dataset after log2 transformation is ", sum(TranscribeGSEA < 0, na.rm=TRUE))
path("# The amount of positive values in TranscribeGSEA dataset after log2 transformation is ", sum(TranscribeGSEA > 0, na.rm=TRUE))
path("# The Normalized TranscribeGSEA dataset has ", sum(TranscribeGSEA==max(TranscribeGSEA), na.rm=TRUE), " Max values and the Max value is ", max(TranscribeGSEA, na.rm=TRUE))
path("# The Normalized TranscribeGSEA dataset has ", sum(TranscribeGSEA==min(TranscribeGSEA), na.rm=TRUE), " min values and the min value is ", min(TranscribeGSEA, na.rm=TRUE))
path("###############################################################################")
path("###############################################################################")
path("# Saving the normalized datasets (by log2 transformation) after low-level processing now ")
print(Sys.time())
d <- Sys.Date()
write.csv(TranscribeGSEA, file=paste("./TranscribeGSEA_Log2Transformed", d, ".csv", sep=""), row.names=TRUE)
LV_Log2 <- TranscribeGSEA
print(head(LV_Log2))
print(head(LV_PEER10_Fudge20_Log2))
print(head(TranscribeCuffdiff))

path("# Making the varied intersections between eGeneLists and expression datasets now ")
LV_Log2_iBaseline <- LV_Log2[intersect(rownames(LV_Log2), imputedBaseline_eGenes), ]
LV_Log2_iIschemia <- LV_Log2[intersect(rownames(LV_Log2), imputedIschemia_eGenes), ]
LV_Log2_Delta     <- LV_Log2[intersect(rownames(LV_Log2), Delta_eGenes), ]
LV_Log2_Ratio     <- LV_Log2[intersect(rownames(LV_Log2), Ratio_eGenes), ]
LV_PEER10_Fudge20_Log2_iBaseline <- LV_PEER10_Fudge20_Log2[intersect(rownames(LV_PEER10_Fudge20_Log2), imputedBaseline_eGenes), ]
LV_PEER10_Fudge20_Log2_iIschemia <- LV_PEER10_Fudge20_Log2[intersect(rownames(LV_PEER10_Fudge20_Log2), imputedIschemia_eGenes), ]
LV_PEER10_Fudge20_Log2_Delta     <- LV_PEER10_Fudge20_Log2[intersect(rownames(LV_PEER10_Fudge20_Log2), Delta_eGenes), ]
LV_PEER10_Fudge20_Log2_Ratio     <- LV_PEER10_Fudge20_Log2[intersect(rownames(LV_PEER10_Fudge20_Log2), Ratio_eGenes), ]
TranscribeCuffdiff_iBaseline <- TranscribeCuffdiff[intersect(rownames(TranscribeCuffdiff), imputedBaseline_eGenes), ]
TranscribeCuffdiff_iIschemia <- TranscribeCuffdiff[intersect(rownames(TranscribeCuffdiff), imputedIschemia_eGenes), ]
TranscribeCuffdiff_Delta     <- TranscribeCuffdiff[intersect(rownames(TranscribeCuffdiff), Delta_eGenes), ]
TranscribeCuffdiff_Ratio     <- TranscribeCuffdiff[intersect(rownames(TranscribeCuffdiff), Ratio_eGenes), ]
LV_Log2_iBaselineDE <- LV_Log2[intersect(rownames(TranscribeCuffdiffSignificant), imputedBaseline_eGenes), ] # 391 eGene/DE
LV_Log2_iIschemiaDE <- LV_Log2[intersect(rownames(TranscribeCuffdiffSignificant), imputedIschemia_eGenes), ] # 325 eGene/DE
LV_Log2_DeltaDE     <- LV_Log2[intersect(rownames(TranscribeCuffdiffSignificant), Delta_eGenes), ]
LV_Log2_RatioDE     <- LV_Log2[intersect(rownames(TranscribeCuffdiffSignificant), Ratio_eGenes), ]
TranscribeCuffdiff_iBaselineDE <- TranscribeCuffdiffSignificant[intersect(rownames(TranscribeCuffdiffSignificant), imputedBaseline_eGenes), ]
TranscribeCuffdiff_iIschemiaDE <- TranscribeCuffdiffSignificant[intersect(rownames(TranscribeCuffdiffSignificant), imputedIschemia_eGenes), ]
TranscribeCuffdiff_DeltaDE     <- TranscribeCuffdiffSignificant[intersect(rownames(TranscribeCuffdiffSignificant), Delta_eGenes), ]
TranscribeCuffdiff_RatioDE     <- TranscribeCuffdiffSignificant[intersect(rownames(TranscribeCuffdiffSignificant), Ratio_eGenes), ]
write.csv(TranscribeCuffdiff_iBaselineDE, file=paste("./TranscribeCuffdiff_iBaselineDE_", d, ".csv", sep=""), row.names=TRUE)
write.csv(TranscribeCuffdiff_iIschemiaDE, file=paste("./TranscribeCuffdiff_iIschemiaDE_", d, ".csv", sep=""), row.names=TRUE)
write.csv(TranscribeCuffdiff_DeltaDE, file=paste("./TranscribeCuffdiff_DeltaDE_", d, ".csv", sep=""), row.names=TRUE)
write.csv(TranscribeCuffdiff_RatioDE, file=paste("./TranscribeCuffdiff_RatioDE_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_Log2_iBaseline, file=paste("./LV_Log2_iBaseline_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_Log2_iIschemia, file=paste("./LV_Log2_iIschemia_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_Log2_Delta, file=paste("./LV_Log2_Delta_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_Log2_Ratio, file=paste("./LV_Log2_Ratio_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_PEER10_Fudge20_Log2_iBaseline, file=paste("./LV_PEER10_Fudge20_Log2_iBaseline_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_PEER10_Fudge20_Log2_iIschemia, file=paste("./LV_PEER10_Fudge20_Log2_iIschemia_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_PEER10_Fudge20_Log2_Delta, file=paste("./LV_PEER10_Fudge20_Log2_Delta_", d, ".csv", sep=""), row.names=TRUE)
write.csv(LV_PEER10_Fudge20_Log2_Ratio, file=paste("./LV_PEER10_Fudge20_Log2_Ratio_", d, ".csv", sep=""), row.names=TRUE)
write.csv(TranscribeCuffdiff_iBaseline, file=paste("./TranscribeCuffdiff_iBaseline_", d, ".csv", sep=""), row.names=TRUE)
write.csv(TranscribeCuffdiff_iIschemia, file=paste("./TranscribeCuffdiff_iIschemia_", d, ".csv", sep=""), row.names=TRUE)
write.csv(TranscribeCuffdiff_Delta, file=paste("./TranscribeCuffdiff_Delta_", d, ".csv", sep=""), row.names=TRUE)
write.csv(TranscribeCuffdiff_Ratio, file=paste("./TranscribeCuffdiff_Ratio_", d, ".csv", sep=""), row.names=TRUE)
colnames(LV_Log2_Ratio)[135] # "B0164V_Pre"
colnames(LV_Log2_Ratio)[136] # "B0002V_Post"
colnames(LV_Log2_Ratio)[268] # "B0164V_Post"
colnames(LV_Log2_Ratio)[269] # NA
colnames(LV_PEER10_Fudge20_Log2_Ratio)[135] # "B0164V_Pre"
colnames(LV_PEER10_Fudge20_Log2_Ratio)[136] # "B0002V_Post"
colnames(LV_PEER10_Fudge20_Log2_Ratio)[268] # NA
colnames(LV_PEER10_Fudge20_Log2_Ratio)[267] # "B0164V_Post"

max(LV_Log2_iBaselineDE) # 13.99293
min(LV_Log2_iBaselineDE) # 0.1375035
LV_Log2_iBaseline_BaselineDE <- LV_Log2_iBaselineDE[,   1:135]
LV_Log2_iBaseline_IschemiaDE <- LV_Log2_iBaselineDE[, 136:268]
LV_Log2_iBaseline_BaselineDE[, 136] <- 15
colnames(LV_Log2_iBaseline_BaselineDE)[136] <- "Boundary"
LV_Log2_iBaselineDEMatrix <- data.matrix(cbind(LV_Log2_iBaseline_BaselineDE, LV_Log2_iBaseline_IschemiaDE))
pdf(file=paste("./LV_Log2_iBaselineDEMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 10, height = 10)
pheatmap(LV_Log2_iBaselineDEMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of 391 Baseline-eGenes/DE using original RNA-Seq (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=0.8, fontsize_row=0.6, color=colorRampPalette(c("white", "blue", "red"))(4096))
dev.off()
max(LV_Log2_iIschemiaDE) # 12.88764
min(LV_Log2_iIschemiaDE) # 0.1375035
LV_Log2_iIschemia_BaselineDE <- LV_Log2_iIschemiaDE[,   1:135]
LV_Log2_iIschemia_IschemiaDE <- LV_Log2_iIschemiaDE[, 136:268]
LV_Log2_iIschemia_BaselineDE[, 136] <- 14
colnames(LV_Log2_iIschemia_BaselineDE)[136] <- "Boundary"
LV_Log2_iIschemiaDEMatrix <- data.matrix(cbind(LV_Log2_iIschemia_BaselineDE, LV_Log2_iIschemia_IschemiaDE))
pdf(file=paste("./LV_Log2_iIschemiaDEMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 10, height = 10)
pheatmap(LV_Log2_iIschemiaDEMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of 325 Ischemia-eGenes/DE using original RNA-Seq (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=0.8, fontsize_row=0.6, color=colorRampPalette(c("white", "firebrick1", "firebrick3"))(4096))
dev.off()

max(LV_Log2_iBaseline) # 16.57207
min(LV_Log2_iBaseline) # 0.1375035
LV_Log2_iBaseline_Baseline <- LV_Log2_iBaseline[,   1:135]
LV_Log2_iBaseline_Ischemia <- LV_Log2_iBaseline[, 136:268]
LV_Log2_iBaseline_Baseline[, 136] <- 17
colnames(LV_Log2_iBaseline_Baseline)[136] <- "Boundary"
LV_Log2_iBaselineMatrix <- data.matrix(cbind(LV_Log2_iBaseline_Baseline, LV_Log2_iBaseline_Ischemia))
pdf(file=paste("./LV_Log2_iBaselineMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 12, height = 15)
pheatmap(LV_Log2_iBaselineMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of Baseline 6140 eGenes using original RNA-Seq (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=0.8, fontsize_row=0.6, color=colorRampPalette(c("white", "blue", "red"))(4096))
dev.off()
max(LV_PEER10_Fudge20_Log2_iBaseline) # 5.197177
min(LV_PEER10_Fudge20_Log2_iBaseline) # 3.928479
LV_PEER10_Fudge20_Log2_iBaseline_Baseline <- LV_PEER10_Fudge20_Log2_iBaseline[,   1:135]
LV_PEER10_Fudge20_Log2_iBaseline_Ischemia <- LV_PEER10_Fudge20_Log2_iBaseline[, 136:267]
LV_PEER10_Fudge20_Log2_iBaseline_Baseline[, 136] <- 6
colnames(LV_PEER10_Fudge20_Log2_iBaseline_Baseline)[136] <- "Boundary"
LV_PEER10_Fudge20_Log2_iBaselineMatrix <- data.matrix(cbind(LV_PEER10_Fudge20_Log2_iBaseline_Baseline, LV_PEER10_Fudge20_Log2_iBaseline_Ischemia))
pdf(file=paste("./LV_PEER10_Fudge20_Log2_iBaselineMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 12, height = 15)
pheatmap(LV_PEER10_Fudge20_Log2_iBaselineMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of Baseline 6140 eGenes using PEER-10 residuals (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=0.8, fontsize_row=0.6, color=colorRampPalette(c("white", "blue", "red"))(4096))
dev.off()

max(LV_Log2_iIschemia) # 14.38324
min(LV_Log2_iIschemia) # 0.1375035
LV_Log2_iIschemia_Baseline <- LV_Log2_iIschemia[,   1:135]
LV_Log2_iIschemia_Ischemia <- LV_Log2_iIschemia[, 136:268]
LV_Log2_iIschemia_Baseline[, 136] <- 15
colnames(LV_Log2_iIschemia_Baseline)[136] <- "Boundary"
LV_Log2_iIschemiaMatrix <- data.matrix(cbind(LV_Log2_iIschemia_Baseline, LV_Log2_iIschemia_Ischemia))
pdf(file=paste("./LV_Log2_iIschemiaMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 12, height = 15)
pheatmap(LV_Log2_iIschemiaMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of Ischemia 5490 eGenes using original RNA-Seq (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=0.8, fontsize_row=0.6, color=colorRampPalette(c("white", "firebrick1", "firebrick3"))(4096))
dev.off()
max(LV_PEER10_Fudge20_Log2_iIschemia) # 5.016221
min(LV_PEER10_Fudge20_Log2_iIschemia) # 4.001477
LV_PEER10_Fudge20_Log2_iIschemia_Baseline <- LV_PEER10_Fudge20_Log2_iIschemia[,   1:135]
LV_PEER10_Fudge20_Log2_iIschemia_Ischemia <- LV_PEER10_Fudge20_Log2_iIschemia[, 136:267]
LV_PEER10_Fudge20_Log2_iIschemia_Baseline[, 136] <- 6
colnames(LV_PEER10_Fudge20_Log2_iIschemia_Baseline)[136] <- "Boundary"
LV_PEER10_Fudge20_Log2_iIschemiaMatrix <- data.matrix(cbind(LV_PEER10_Fudge20_Log2_iIschemia_Baseline, LV_PEER10_Fudge20_Log2_iIschemia_Ischemia))
pdf(file=paste("./LV_PEER10_Fudge20_Log2_iIschemiaMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 12, height = 15)
pheatmap(LV_PEER10_Fudge20_Log2_iIschemiaMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of Ischemia 5490 eGenes using PEER-10 residuals (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=0.8, fontsize_row=0.6, color=colorRampPalette(c("white", "firebrick1", "firebrick3"))(4096))
dev.off()
















max(LV_Log2_Delta) # 12.28021
min(LV_Log2_Delta) # 0.1375035
LV_Log2_Delta_Baseline <- LV_Log2_Delta[,   1:135]
LV_Log2_Delta_Ischemia <- LV_Log2_Delta[, 136:268]
LV_Log2_Delta_Baseline[, 136] <- 13
colnames(LV_Log2_Delta_Baseline)[136] <- "Boundary"
LV_Log2_DeltaMatrix <- data.matrix(cbind(LV_Log2_Delta_Baseline, LV_Log2_Delta_Ischemia))
pdf(file=paste("./LV_Log2_DeltaMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 10, height = 8)
pheatmap(LV_Log2_DeltaMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of 32 Response_Delta eGenes using original RNA-Seq (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=1, color=colorRampPalette(c("white","yellow","orange","firebrick2","red4"))(4096))
dev.off()
max(LV_PEER10_Fudge20_Log2_Delta) # 5.001907
min(LV_PEER10_Fudge20_Log2_Delta) # 4.183503
LV_PEER10_Fudge20_Log2_Delta_Baseline <- LV_PEER10_Fudge20_Log2_Delta[,   1:135]
LV_PEER10_Fudge20_Log2_Delta_Ischemia <- LV_PEER10_Fudge20_Log2_Delta[, 136:267]
LV_PEER10_Fudge20_Log2_Delta_Baseline[, 136] <- 6
colnames(LV_PEER10_Fudge20_Log2_Delta_Baseline)[136] <- "Boundary"
LV_PEER10_Fudge20_Log2_DeltaMatrix <- data.matrix(cbind(LV_PEER10_Fudge20_Log2_Delta_Baseline, LV_PEER10_Fudge20_Log2_Delta_Ischemia))
pdf(file=paste("./LV_PEER10_Fudge20_Log2_DeltaMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 10, height = 8)
pheatmap(LV_PEER10_Fudge20_Log2_DeltaMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of 32 Response_Delta eGenes using PEER-10 residuals (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=1, color=colorRampPalette(c("white", "blue", "red"))(4096))
dev.off()

max(LV_Log2_Ratio) # 8.483461
LV_Log2_Ratio_Baseline <- LV_Log2_Ratio[,   1:135]
LV_Log2_Ratio_Ischemia <- LV_Log2_Ratio[, 136:268]
LV_Log2_Ratio_Baseline[, 136] <- 9
colnames(LV_Log2_Ratio_Baseline)[136] <- "Boundary"
LV_Log2_RatioMatrix <- data.matrix(cbind(LV_Log2_Ratio_Baseline, LV_Log2_Ratio_Ischemia))
pdf(file=paste("./LV_Log2_RatioMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 10, height = 8)
pheatmap(LV_Log2_RatioMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of 19 Response_Ratio eGenes using original RNA-Seq (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=1, color=colorRampPalette(c("white","darkorange","orangered","red2"))(4096))
dev.off()
max(LV_PEER10_Fudge20_Log2_Ratio) # 4.820488
min(LV_PEER10_Fudge20_Log2_Ratio) # 4.277816
LV_PEER10_Fudge20_Log2_Ratio_Baseline <- LV_PEER10_Fudge20_Log2_Ratio[,   1:135]
LV_PEER10_Fudge20_Log2_Ratio_Ischemia <- LV_PEER10_Fudge20_Log2_Ratio[, 136:267]
LV_PEER10_Fudge20_Log2_Ratio_Baseline[, 136] <- 5
colnames(LV_PEER10_Fudge20_Log2_Ratio_Baseline)[136] <- "Boundary"
LV_PEER10_Fudge20_Log2_RatioMatrix <- data.matrix(cbind(LV_PEER10_Fudge20_Log2_Ratio_Baseline, LV_PEER10_Fudge20_Log2_Ratio_Ischemia))
pdf(file=paste("./LV_PEER10_Fudge20_Log2_RatioMatrix_eGenesHeatMap.pdf", sep=""), paper="special", width = 10, height = 8)
pheatmap(LV_PEER10_Fudge20_Log2_RatioMatrix, cluster_rows=TRUE, cluster_cols=FALSE, 
 main=paste("Expression Heatmap of Response_Ratio eGenes after PEER10 regressOut (left-Baseline right-Ischemia)", sep=""), 
 fontsize_col=1, color=colorRampPalette(c("white", "blue", "red"))(4096))
dev.off()





















 
iDmatrix <- data.matrix(iDFrame)
  pdf(file=paste("./imputedDelta_residualPEER", knum, "_eGenesHeatMap.pdf", sep=""), paper="special", width = 9, height = 7)
   pheatmap(iDmatrix, legend=TRUE, legend_breaks=c(0,1), legend_labels=c("Insignificant","Significant"), 
   cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=FALSE, cellwidth=36, cellheight=5, 
   main=paste("eGenes Heatmap of imputedDelta_residualPEER-", knum, " under 2nd Pass of PEER-K", sep=""), 
   color=colorRampPalette(c("navy","white","firebrick3"))(4096), fontsize=10, fontsize_row=6, fontsize_col=10, fontsize_number=8)
  dev.off()




# length(levels(factor(Delta_Response$probeid)))
          Delta_eGenes <- levels(factor(Delta_Response$probeid))
          Ratio_eGenes <- levels(factor(Ratio_Response$probeid))
imputedBaseline_eGenes <- levels(factor(imputedBaseline_FDR.05.gQTLstats_PEER10$probeid))
imputedIschemia_eGenes <- levels(factor(imputedIschemia_FDR.05.gQTLstats_PEER10$probeid))
save(Delta_eGenes, Ratio_eGenes, imputedBaseline_eGenes, imputedIschemia_eGenes, file="./eGeneList.rda")
eGeneList <- list(Delta_eGene=Delta_eGenes, Ratio_eGene=Ratio_eGenes, Baseline_eGene=imputedBaseline_eGenes, Ischemia_eGene=imputedIschemia_eGenes)
names(eGeneList)
# test.df <- do.call("rbind", lapply(eGeneList, as.data.frame)) 
# writeLines(unlist(eGeneList), "newfile.txt")
write.table(eGeneList$Delta_eGene,    file="./Delta_eGene.txt",    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(eGeneList$Ratio_eGene,    file="./Ratio_eGene.txt",    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(eGeneList$Baseline_eGene, file="./Baseline_eGene.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(eGeneList$Ischemia_eGene, file="./Ischemia_eGene.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

library("VennDiagram")
library("scatterplot3d")

# Reference four-set diagram
venn.plot <- draw.quad.venn(
area1 = length(imputedBaseline_eGenes),
area2 = length(imputedIschemia_eGenes),
area3 = length(Delta_eGenes),
area4 = length(Ratio_eGenes),
n12   = length(intersect(imputedBaseline_eGenes,imputedIschemia_eGenes)),
n13   = length(intersect(imputedBaseline_eGenes,Delta_eGenes)),
n14   = length(intersect(imputedBaseline_eGenes,Ratio_eGenes)),
n23   = length(intersect(imputedIschemia_eGenes,Delta_eGenes)),
n24   = length(intersect(imputedIschemia_eGenes,Ratio_eGenes)),
n34   = length(intersect(Delta_eGenes,Ratio_eGenes)),
n123  = length(intersect(intersect(imputedBaseline_eGenes,imputedIschemia_eGenes),Delta_eGenes)),
n124  = length(intersect(intersect(imputedBaseline_eGenes,imputedIschemia_eGenes),Ratio_eGenes)),
n134  = length(intersect(intersect(imputedBaseline_eGenes,Delta_eGenes),Ratio_eGenes)),
n234  = length(intersect(intersect(imputedIschemia_eGenes,Delta_eGenes),Ratio_eGenes)),
n1234 = length(intersect(intersect(imputedBaseline_eGenes,imputedIschemia_eGenes),intersect(Delta_eGenes,Ratio_eGenes))),
category = c("Baseline(6140)", "Ischemia(5490)", "Delta(32)", "Ratio(19)"),
fill = c("orange", "red", "green", "blue"),
lty = "dashed",
cex = 3,
cat.cex = 1.6,
cat.col = c("orange", "red", "green", "blue") );
grid.draw(venn.plot);
grid.newpage();
# tiff(filename = "Quad_Venn_diagram.tiff", compression = "lzw");
pdf(file="./Quad_VennDiagram.pdf", paper="special", width = 12, height = 8)
 grid.draw(venn.plot);
dev.off();

# Reference three-set diagram
Response_eGenes <- union(Delta_eGenes, Ratio_eGenes)
venn.plot <- draw.triple.venn(
area1 = length(imputedBaseline_eGenes),
area2 = length(imputedIschemia_eGenes),
area3 = length(Response_eGenes),
n12   = length(intersect(imputedBaseline_eGenes,imputedIschemia_eGenes)),
n23   = length(intersect(imputedIschemia_eGenes,Response_eGenes)),
n13   = length(intersect(imputedBaseline_eGenes,Response_eGenes)),
n123  = length(intersect(intersect(imputedBaseline_eGenes,imputedIschemia_eGenes),Response_eGenes)),
category = c("Baseline(6140)", "Ischemia(5490)", "Response(51)"),
fill     = c("red", "blue", "green"),
lty      = "blank",
cex      = 4,   # The font size of the number in the circle
cat.cex  = 1.6, # The font size of Baseline_F ....etc
cat.col  = c("red", "blue", "green") );
grid.draw(venn.plot);
grid.newpage();
# tiff(filename = "Triple_Venn_diagram.tiff", compression = "lzw");
pdf(file="./Triple_VennDiagram.pdf", paper="special", width = 12, height = 8)
 grid.draw(venn.plot);
dev.off();





Delta_Response <- Delta_Response[, -1]
Ratio_Response <- Ratio_Response[, -1]
names(HumanLV_Exp) <- c("Gene_ID","RPKM_Count")
Delta_HumanLV    <- intersect(levels(factor(Delta_Response$probeid)), levels(factor(HumanLV_Exp$Gene_ID)))
Delta_MouseHeart <- intersect(levels(factor(Delta_Response$probeid)), levels(factor(MouseHeart_Exp$Gene)))
# Delta_Response[Delta_Response$probeid=="COL4A1", ]

HumanLV_Exp2 <- HumanLV_Exp[HumanLV_Exp$Gene_ID %in% Delta_HumanLV, ]
chisqG <- c(1:length(Delta_HumanLV))
pifdr  <- c(1:length(Delta_HumanLV))
rN <- length(Delta_HumanLV)
for(i in seq(rN)) {
	 test <- Delta_Response[Delta_Response$probeid %in% Delta_HumanLV[i], ]
   chisqG[i] <- max(test$chisq)
   pifdr[i]  <- test$fdr[min(which(test$chisq==chisqG[i]))]
}
HumanLV_Exp2$chisq <- chisqG
HumanLV_Exp2$pifdr <- pifdr
Delta_HumanLV <- HumanLV_Exp2
rm(HumanLV_Exp2)
Delta_HumanLV <- Delta_HumanLV[order(Delta_HumanLV$pifdr), ]
Delta_HumanLV <- Delta_HumanLV[order(Delta_HumanLV$chisq, decreasing = TRUE), ]
write.csv(Delta_HumanLV, file="./Delta_HumanLV.csv", row.names=FALSE)
library("scatterplot3d")
pdf(file="./deltaResponseHumanLV.pdf", paper="special", width = 9, height = 7)
with(Delta_HumanLV, {
    s3d <- scatterplot3d(chisq, pifdr, RPKM_Count,                   # x y and z axis
                       main = "29 deltaResponse eGenes after filtered by Human LV",
                       highlight.3d = TRUE,
                       grid = TRUE,
                       type = "h",
                       col.axis = "blue",
                       col.grid = "lightblue",
                       angle = 65,
                       scale.y = 0.6, 
                       pch = 20,
                       xlab="ChiSQ",
                       ylab="piFDR (plug-in FDR)",
                       zlab="Expression Count (RPKM)")
    s3d.coords <- s3d$xyz.convert(chisq, pifdr, RPKM_Count) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,                        # x and y coordinates                       
         #labels=row.names(mtcars),                         # text to plot
          labels=Gene_ID,
          cex=.5, pos=4)                                    # shrink text 50% and place to right of points)
})
dev.off()

MouseHeart_Exp2 <- MouseHeart_Exp[MouseHeart_Exp$Gene %in% Delta_MouseHeart, ]
MouseHeart_Exp2 <- MouseHeart_Exp2[, -c(2,4)]
names(MouseHeart_Exp2) <- c("Gene_ID","mHeartPercentile")
chisqG <- c(1:length(Delta_MouseHeart))
pifdr  <- c(1:length(Delta_MouseHeart))
rN <- length(Delta_MouseHeart)
for(i in seq(rN)) {
	 test <- Delta_Response[Delta_Response$probeid %in% Delta_MouseHeart[i], ]
   chisqG[i] <- max(test$chisq)
   pifdr[i]  <- test$fdr[min(which(test$chisq==chisqG[i]))]
}
MouseHeart_Exp2$chisq <- chisqG
MouseHeart_Exp2$pifdr <- pifdr
Delta_MouseHeart <- MouseHeart_Exp2
rm(MouseHeart_Exp2)


Delta_MouseHeart <- Delta_MouseHeart[order(Delta_MouseHeart$pifdr), ]
Delta_MouseHeart <- Delta_MouseHeart[order(Delta_MouseHeart$chisq, decreasing = TRUE), ]
write.csv(Delta_MouseHeart, file="./Delta_MouseHeart.csv", row.names=FALSE)
library("scatterplot3d")
pdf(file="./deltaResponseMouseHeart.pdf", paper="special", width = 9, height = 7)
with(Delta_MouseHeart, {
    s3d <- scatterplot3d(chisq, pifdr, mHeartPercentile,                   # x y and z axis
                       main = "3 deltaResponse eGenes after filtered by Mouse top75%-tile expression",
                       highlight.3d = TRUE,
                       grid = TRUE,
                       type = "h",
                       col.axis = "blue",
                       col.grid = "lightblue",
                       angle = 55,
                       scale.y = 0.6, 
                       pch = 20,
                       xlab="ChiSQ",
                       ylab="piFDR (plug-in FDR)",
                       zlab="Expression Percentile Rank (>75%)")
    s3d.coords <- s3d$xyz.convert(chisq, pifdr, mHeartPercentile) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,                        # x and y coordinates                       
         #labels=row.names(mtcars),                         # text to plot
          labels=Gene_ID,
          cex=.5, pos=4)                                    # shrink text 50% and place to right of points)
})
dev.off()

response_intersect <- intersect(levels(factor(Delta_Response$probeid)), levels(factor(Ratio_Response$probeid)))

















Ratio_HumanLV    <- intersect(levels(factor(Ratio_Response$probeid)), levels(factor(HumanLV_Exp$Gene_ID)))
Ratio_MouseHeart <- intersect(levels(factor(Ratio_Response$probeid)), levels(factor(MouseHeart_Exp$Gene)))
# Ratio_Response[Ratio_Response$probeid=="COL4A1", ]

HumanLV_Exp2 <- HumanLV_Exp[HumanLV_Exp$Gene_ID %in% Ratio_HumanLV, ]
chisqG <- c(1:length(Ratio_HumanLV))
pifdr  <- c(1:length(Ratio_HumanLV))
rN <- length(Ratio_HumanLV)
for(i in seq(rN)) {
	 test <- Ratio_Response[Ratio_Response$probeid %in% Ratio_HumanLV[i], ]
   chisqG[i] <- max(test$chisq)
   pifdr[i]  <- test$fdr[min(which(test$chisq==chisqG[i]))]
}
HumanLV_Exp2$chisq <- chisqG
HumanLV_Exp2$pifdr <- pifdr
Ratio_HumanLV <- HumanLV_Exp2
rm(HumanLV_Exp2)
Ratio_HumanLV <- Ratio_HumanLV[order(Ratio_HumanLV$pifdr), ]
Ratio_HumanLV <- Ratio_HumanLV[order(Ratio_HumanLV$chisq, decreasing = TRUE), ]
write.csv(Ratio_HumanLV, file="./Ratio_HumanLV.csv", row.names=FALSE)
library("scatterplot3d")
pdf(file="./ratioResponseHumanLV.pdf", paper="special", width = 9, height = 7)
with(Ratio_HumanLV, {
    s3d <- scatterplot3d(chisq, pifdr, RPKM_Count,                   # x y and z axis
                       main = "18 ratioResponse eGenes after filtered by Human LV",
                       highlight.3d = TRUE,
                       grid = TRUE,
                       type = "h",
                       col.axis = "blue",
                       col.grid = "lightblue",
                       angle = 65,
                       scale.y = 0.6, 
                       pch = 20,
                       xlab="ChiSQ",
                       ylab="piFDR (plug-in FDR)",
                       zlab="Expression Count (RPKM)")
    s3d.coords <- s3d$xyz.convert(chisq, pifdr, RPKM_Count) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,                        # x and y coordinates                       
         #labels=row.names(mtcars),                         # text to plot
          labels=Gene_ID,
          cex=.5, pos=4)                                    # shrink text 50% and place to right of points)
})
dev.off()

MouseHeart_Exp2 <- MouseHeart_Exp[MouseHeart_Exp$Gene %in% Ratio_MouseHeart, ]
MouseHeart_Exp2 <- MouseHeart_Exp2[, -c(2,4)]
names(MouseHeart_Exp2) <- c("Gene_ID","mHeartPercentile")
chisqG <- c(1:length(Ratio_MouseHeart))
pifdr  <- c(1:length(Ratio_MouseHeart))
rN <- length(Ratio_MouseHeart)
for(i in seq(rN)) {
	 test <- Ratio_Response[Ratio_Response$probeid %in% Ratio_MouseHeart[i], ]
   chisqG[i] <- max(test$chisq)
   pifdr[i]  <- test$fdr[min(which(test$chisq==chisqG[i]))]
}
MouseHeart_Exp2$chisq <- chisqG
MouseHeart_Exp2$pifdr <- pifdr
Ratio_MouseHeart <- MouseHeart_Exp2
rm(MouseHeart_Exp2)


Ratio_MouseHeart <- Ratio_MouseHeart[order(Ratio_MouseHeart$pifdr), ]
Ratio_MouseHeart <- Ratio_MouseHeart[order(Ratio_MouseHeart$chisq, decreasing = TRUE), ]
write.csv(Ratio_MouseHeart, file="./Ratio_MouseHeart.csv", row.names=FALSE)
library("scatterplot3d")
pdf(file="./ratioResponseMouseHeart.pdf", paper="special", width = 9, height = 7)
with(Ratio_MouseHeart, {
    s3d <- scatterplot3d(chisq, pifdr, mHeartPercentile,                   # x y and z axis
                       main = "6 ratioResponse eGenes after filtered by Mouse top75%-tile expression",
                       highlight.3d = TRUE,
                       grid = TRUE,
                       type = "h",
                       col.axis = "blue",
                       col.grid = "lightblue",
                       angle = 55,
                       scale.y = 0.6, 
                       pch = 20,
                       xlab="ChiSQ",
                       ylab="piFDR (plug-in FDR)",
                       zlab="Expression Percentile Rank (>75%)")
    s3d.coords <- s3d$xyz.convert(chisq, pifdr, mHeartPercentile) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,                        # x and y coordinates                       
         #labels=row.names(mtcars),                         # text to plot
          labels=Gene_ID,
          cex=.5, pos=4)                                    # shrink text 50% and place to right of points)
})
dev.off()



















