path Aug 29 11:19:32 EDT 2016
Intersecting the compath shared SNPs by their SNPs loci now
path Aug 29 11:20:44 EDT 2016
The intersected SNPs between two files contains 111774 rows
path Aug 29 11:20:44 EDT 2016
The intersected SNPs between two files contains 2513 columns
path Aug 29 11:20:44 EDT 2016
Finished extracting the compath shared SNPs between two files now
======================================================================================================================
setwd("D:/BWHMS/Meetings02/2016-09-01")
library(data.table)
commShared_chr1 <- fread("commShared_chr1.vcf")
 rsSNP_trans142 <- fread("rsSNPs_Transcribe142_Chr01_nonimputed.vcf")
commShared_chr1 <- as.data.frame(commShared_chr1)
 rsSNP_trans142 <- as.data.frame(rsSNP_trans142)
write.table(commShared_chr1$ID, file="./compathSNPs.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
# Error message from previous merging:
# Found 444 SNPs that do not match in terms of allele codes
# Might include strand flips, although flipped A/T and C/G SNPs will be undetected)
# Writing problem SNPs to [ combined.missnp ]
# So we have to remove these 444 SNPs from the compath shared SNPs
compathSNPs <- as.data.frame(fread("compathSNPs.txt", header=FALSE))
removedSNP <- as.data.frame(fread("combined.missnp.txt", header=FALSE))
newCompathSNPs <- as.data.frame(setdiff(compathSNPs$V1, removedSNP$V2))
names(newCompathSNPs) <- c("SNP_ID")
write.table(newCompathSNPs$SNP_ID, file="./newCompathSNPs.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

======================================================================================================================
The following session has to be run under UNIX BASH Shell 
======================================================================================================================
qrsh -l lx6 -l mem_free=92G
cd /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160829PCA

# Convert the VCF back to PED/MAP files
/app/vcftools-0.1.12a@i86-rhel6.0/vcftools_0.1.12a/bin/vcftools --vcf ./rsSNPs_Transcribe142_Chr01_nonimputed.vcf \
--plink --remove-filtered-geno-all --out rsSNPs_Transcribe142_Chr01_nonimputed.recode

# generated the BED/BIM/FAM files
/app/plink-2@i86-rhel6.0/bin/plink --noweb --file rsSNPs_Transcribe142_Chr01_nonimputed.recode --make-bed --out rsSNPs_Transcribe142_Chr01_nonimputed.recode

# extract only the compath SNPs
plink --file ALL.chr1.phase3.20130502 \
      --extract newCompathSNPs.txt \
      --recode --out All_1000_genome
plink --file rsSNPs_Transcribe142_Chr01_nonimputed.recode \
      --extract newCompathSNPs.txt \
      --recode --out rsSNPs_Trans142

# merge the two files
plink --file rsSNPs_Trans142 \
      --merge All_1000_genome.ped \
              All_1000_genome.map \
      --recode --out combined
plink --noweb \
      --file combined \
      --make-bed \
      --out combined
       
cd /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source
for chr in `seq 25`; do cp /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160829PCA/combined.ped /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source/chr${chr}/PCA_chr${chr}.ped; done
for chr in `seq 25`; do cp /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160829PCA/combined.map /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source/chr${chr}/PCA_chr${chr}.map; done
for chr in `seq 25`; do cp /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160829PCA/combined.bed /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source/chr${chr}/PCA_chr${chr}.bed; done
for chr in `seq 25`; do cp /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160829PCA/combined.bim /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source/chr${chr}/PCA_chr${chr}.bim; done
for chr in `seq 25`; do cp /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160829PCA/combined.fam /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source/chr${chr}/PCA_chr${chr}.fam; done

#run PCA#
# qrsh -l lx6 -l mem_free=92G
source /proj/relibs/relib00/python/2.7.3_x86_64_CentOS6.5/bin/activate
cd /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Results
cdnm_pca_pipeline --parallel=100 /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source/
cdnm_pca_pipeline --out-base=Trans1000_20160829_combined /udd/retwc/TRANSCRIBE/analyses/retwc/20160829PCA/PCA_Source/
make
# Finished running PCA#




library(Matrix)
library(irlba)
library(plyr)
library(threejs)
indEthnic = readLines(pipe("zcat ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | sed -n /^#CHROM/p | tr '\t' '\n' | tail -n +10"))
save(indEthnic, file="./indEthnic.rda")
load("indEthnic.rda")
ped = read.table("./20130606_g1k.ped",sep="\t",header=TRUE,row.names=2)[indEthnic,6,drop=FALSE] # Without "drop=FALSE" the ped won't become data frame
pop = read.table("./20131219.populations.tsv",sep="\t",header=TRUE)
pop = pop[1:26,]
super = pop[,3]
names(super) = pop[,2]
super = factor(super)

PCs20 <- read.table("D:/BWHMS/Meetings02/2016-08-25/20160825PCA/PCA_Results/Trans1000_20160825_combined.ProPC.coord.txt", header = TRUE, stringsAsFactors = FALSE)
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")

scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.5)

newPC20 <- merge(PCs20, ped, by="popID", all.x=TRUE)

===========================================================================================================================================================
New PCA (Transcribe120+100GP32504) results  (2016-08-31)
===========================================================================================================================================================
library(Matrix)
library(data.table)
library(irlba)
library(plyr)
library(threejs)

newPC20 <- as.data.frame(fread("./PC20_Transcribe142_SNP50K.txt"))
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")
scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.5)
test <- subset(newPC20, newPC20$col!="yellow", select=-c(Superpopulation))
test <- subset(test, test$col!="dimgrey", select=-c(Population))
newPC20 <- test
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")
scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.5)
test <- subset(newPC20, newPC20$col!="sienna")
newPC20 <- test
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")
scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.3)
test <- subset(newPC20, newPC20$col!="cyan")
newPC20 <- test
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")
scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.6)

newPC20 <- as.data.frame(fread("./PC20_Transcribe142_SNP50K.txt"))
test <- subset(newPC20, newPC20$col!="yellow")
test <- subset(test, test$col!="dimgrey")
test <- subset(test, test$col!="sienna")
test <- subset(test, test$col!="cyan")
levels(factor(test$Population))
test[which(test$Population=="FIN"),8] <- "gold"
test[which(test$Population=="GBR"),8] <- "greenyellow"
test[which(test$Population=="IBS"),8] <- "lightsalpath"
newPC20 <- test
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")
scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.6)

PCs20 <- read.table("D:/BWHMS/Meetings02/2016-09-01/PCA_Results/Trans1000_20160831_combined120.ProPC.coord", header = TRUE, stringsAsFactors = FALSE)


===========================================================================================================================================================
New PCA (Transcribe120+100GP32504) results  (2016-09-08)
===========================================================================================================================================================
library(Matrix)
library(data.table)
library(irlba)
library(plyr)
library(threejs)

PCs20 <- read.table("D:/BWHMS/Meetings02/2016-09-08/PCA_111329SNPs+2646Subjects/Trans1000_20160825_combined.ProPC.coord", header = TRUE, stringsAsFactors = FALSE)
load("D:/BWHMS/Meetings02/2016-09-01/1000GP3_Ethnics.rda")
PCs20 <- PCs20[, -c(2:5, 9:25)]
newPC20 <- merge(PCs20, ped, by="popID", all.x=TRUE)
newPC20[1:142, 7] <- "blue"
newPC20[47, 7] <- "red"
newPC20[72, 7] <- "red"
newPC20[137, 7] <- "red"
write.table(newPC20, file="./PC20_100K.txt", quote=TRUE, row.names=FALSE, col.names=TRUE, sep="\t")
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")
scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.5)

PCs20 <- read.table("D:/BWHMS/Meetings02/2016-09-01/PCA_Trans120/Transcribe120_AllChr22_imputed_MAF0.15passed.ProPC.coord", header = TRUE, stringsAsFactors = FALSE)
newPC20 <- PCs20[, -c(2:5, 9:25)]
newPC20$col <- "blue"
x <- newPC20$PC1
y <- newPC20$PC2
z <- newPC20$PC3
dataPCA <- cbind(x,y,z) # matrix format
colnames(dataPCA) <- c("PC1","PC2","PC3")
scatterplot3js(dataPCA, color=newPC20$col, flip.y=FALSE, size=0.5)

rR <- as.data.frame(fread("D:/BWHMS/Meetings02/2016-09-08/ResponseQTL/imputed_Ratio/ResponseQTL_Ratio_imputed_ratioResidual_FDR.05_PEER-10.csv"))
dR <- as.data.frame(fread("D:/BWHMS/Meetings02/2016-09-08/ResponseQTL/imputed_Delta/ResponseQTL_Delta_imputed_deltaResidual_FDR.05_PEER-10.csv"))
