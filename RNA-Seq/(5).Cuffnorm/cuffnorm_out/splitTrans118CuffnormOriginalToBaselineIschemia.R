
( awk 'BEGIN {OFS="\t"} {print $4, $7}' ./genes.attr_table ) | ( awk -F\: 'BEGIN {OFS="\t"} {print $1, $2}')  > temp1
 sed '1,1d' ./temp1 > temp2 #remove the first row
( awk 'BEGIN {OFS="\t"} {print $3}' ./temp2 ) | ( sed 's/-/\t/g' ) > temp3 
#Extract the genes coordinates which is the $3 and remove the '-' and save it as temp3 
( awk 'BEGIN {OFS="\t"} {print $1, $2}' ./temp2 ) | ( sed 's/ /\t/g' ) > temp4
( paste temp4 temp3 ) > geneCoordinates
# ( awk -F\; 'BEGIN {OFS="\t"} {print $3, $1}' ./mergedID.gtf ) | ( awk '{print $12, $3, $9}' ) | ( sed -n 's/"//gp' ) | awk '!seen[$0]++' > geneid.mergedID
####################################################################################################################
#   When finished running the Cuffnorm, go into the Cuffnorm output folder named "cuffnorm_out"                    #
#   Run a BASH script named "cuffnormID_transcribeID.HMS.lsf" to convert the Cuffnorm ID to our real patient ID    #
#   The script will generate a file named "transcribe118.RNA-seq.cuffnormOriginal"                                 #
#   Read this file into R environment and start the following procedure to split the Baseline and Ischemia dataset #
####################################################################################################################
# require("gQTLstats")
# require("gQTLBase")
require("VariantAnnotation")
require("GenomicRanges")
require("GenomeInfoDb")
require("Rsamtools")
require(data.table)
require("ggplot2")
require("reshape2")
require("graphics")
require("biomaRt")
require("tools")
listMarts()
base.dir <- getwd()
odir <- getwd()
set.seed(12345)

# sampleID <- as.data.frame(fread("./sampleID.txt", header=FALSE))
# transcribe118 <- as.data.frame(fread("./genes.fpkm_table")) # errors happaned
geneID <- read.table("./geneCoordinates", header=FALSE, stringsAsFactors=FALSE)
mergedID <- read.table("./geneid.mergedID", header=TRUE, stringsAsFactors=FALSE)
names(geneID) <- c("ugene_id","chr","start","end")
mergedID$chr <- NULL; mergedID$start <- NULL; mergedID$end <- NULL;
# the '-' caused many troubles while merging the two datasets
# therefore we have to make unique gene_id first
#
# > test <- merge(geneID, mergedID, by="ugene_id", all.x=TRUE)
# > sum(is.na(test$strand))
# [1] 581
# > sum(is.na(mergedID$strand))
# [1] 0
# > testNA <- test[which(is.na(test$strand)),]
# > head(testNA)
#         ugene_id   chr     start       end strand
# 1    1/2-SBSRNA4  chr4 110351118 110461615   <NA>
# 3       A1BG-AS1 chr19  58858171  58874214   <NA>
# 33      AATK-AS1 chr17  79008946  79156964   <NA>
# 88    ABHD11-AS1  chr7  73149398  73150330   <NA>
# 93  ABHD14A-ACY1  chr3  52009041  52023218   <NA>
# 297  ADAMTS9-AS2  chr3  64501330  64997143   <NA>
# > 
geneID$ugene_id <- make.names(geneID$ugene_id, unique=TRUE)
mergedID$ugene_id <- make.names(mergedID$ugene_id, unique=TRUE)
# > test <- merge(geneID, mergedID, by="ugene_id", all.x=TRUE)
# > sum(is.na(test$strand))
# [1] 0 # It proved that there is no NA value existed
# > 
mergedID <- merge(geneID, mergedID, by="ugene_id", all.x=TRUE)
mergedID[,6] <- mergedID$start; mergedID[,7] <- mergedID$end;
mergedID <- mergedID[, -c(3:4)]
names(mergedID) <- c("ugene_id","chr","strand","start","end")
write.csv(mergedID,  file="./mergedID.csv", row.names=FALSE)

trans118 <- read.table("transcribe118.RNA-seq.cuffnormOriginal", header = TRUE, stringsAsFactors = FALSE) #23079 x 237
trans118$tracking_id <- make.names(trans118$tracking_id, unique=TRUE)
rownames(trans118) <- trans118[, 1]
trans118 <- trans118[, -1]
write.table(trans118, file="./transcribe118.RNA-seq.cuffnorm.tab", quote = TRUE, row.names = TRUE, col.names = NA, sep = "\t")
write.table(trans118, file="./transcribe118.RNA-seq.cuffnorm.tab2", quote = TRUE, row.names = TRUE, col.names = NA, sep = "\t")
##  names(trans118)[names(trans118)=="X"] <- "ugene_id"
##> trans118[1:2, 118:119]
##               B0152V_post B0041V_pre
##  1/2-SBSRNA4    1.196820   0.639935
##  A1BG           0.160783   0.750462
LV118post <- trans118[, 1:118]
LV118pre  <- trans118[, -c(1:118)]
LV118pre  <-  LV118pre[,order(names(LV118pre))]
LV118post <- LV118post[,order(names(LV118post))]
#> names(LV118pre)
#  [1] "B0003V_pre" "B0004V_pre" "B0006V_pre" "B0008V_pre" "B0009V_pre" "B0012V_pre" "B0013V_pre"
#  [8] "B0014V_pre" "B0015V_pre" "B0017V_pre" "B0018V_pre" "B0021V_pre" "B0022V_pre" "B0024V_pre"
# [15] "B0027V_pre" "B0028V_pre" "B0029V_pre" "B0030V_pre" "B0031V_pre" "B0032V_pre" "B0033V_pre"
# [22] "B0034V_pre" "B0035V_pre" "B0037V_pre" "B0039V_pre" "B0041V_pre" "B0042V_pre" "B0043V_pre"
# [29] "B0044V_pre" "B0045V_pre" "B0046V_pre" "B0047V_pre" "B0048V_pre" "B0049V_pre" "B0050V_pre"
# [36] "B0051V_pre" "B0052V_pre" "B0054V_pre" "B0055V_pre" "B0056V_pre" "B0060V_pre" "B0061V_pre"
# [43] "B0062V_pre" "B0063V_pre" "B0064V_pre" "B0065V_pre" "B0066V_pre" "B0067V_pre" "B0068V_pre"
# [50] "B0069V_pre" "B0070V_pre" "B0071V_pre" "B0074V_pre" "B0076V_pre" "B0078V_pre" "B0079V_pre"
# [57] "B0080V_pre" "B0081V_pre" "B0082V_pre" "B0086V_pre" "B0088V_pre" "B0089V_pre" "B0092V_pre"
# [64] "B0093V_pre" "B0096V_pre" "B0098V_pre" "B0099V_pre" "B0100V_pre" "B0102V_pre" "B0103V_pre"
# [71] "B0104V_pre" "B0105V_pre" "B0106V_pre" "B0107V_pre" "B0108V_pre" "B0109V_pre" "B0110V_pre"
# [78] "B0111V_pre" "B0112V_pre" "B0113V_pre" "B0114V_pre" "B0115V_pre" "B0116V_pre" "B0118V_pre"
# [85] "B0120V_pre" "B0122V_pre" "B0123V_pre" "B0124V_pre" "B0125V_pre" "B0126V_pre" "B0127V_pre"
# [92] "B0129V_pre" "B0132V_pre" "B0133V_pre" "B0135V_pre" "B0136V_pre" "B0137V_pre" "B0138V_pre"
# [99] "B0139V_pre" "B0140V_pre" "B0141V_pre" "B0142V_pre" "B0143V_pre" "B0144V_pre" "B0145V_pre"
#[106] "B0149V_pre" "B0150V_pre" "B0151V_pre" "B0152V_pre" "B0153V_pre" "B0154V_pre" "B0156V_pre"
#[113] "B0157V_pre" "B0159V_pre" "B0160V_pre" "B0162V_pre" "B0163V_pre" "B0164V_pre"
#> names(LV118post)
#  [1] "B0003V_post" "B0004V_post" "B0006V_post" "B0008V_post" "B0009V_post" "B0012V_post" "B0013V_post"
#  [8] "B0014V_post" "B0015V_post" "B0017V_post" "B0018V_post" "B0021V_post" "B0022V_post" "B0024V_post"
# [15] "B0027V_post" "B0028V_post" "B0029V_post" "B0030V_post" "B0031V_post" "B0032V_post" "B0033V_post"
# [22] "B0034V_post" "B0035V_post" "B0037V_post" "B0039V_post" "B0041V_post" "B0042V_post" "B0043V_post"
# [29] "B0044V_post" "B0045V_post" "B0046V_post" "B0047V_post" "B0048V_post" "B0049V_post" "B0050V_post"
# [36] "B0051V_post" "B0052V_post" "B0054V_post" "B0055V_post" "B0056V_post" "B0060V_post" "B0061V_post"
# [43] "B0062V_post" "B0063V_post" "B0064V_post" "B0065V_post" "B0066V_post" "B0067V_post" "B0068V_post"
# [50] "B0069V_post" "B0070V_post" "B0071V_post" "B0074V_post" "B0076V_post" "B0078V_post" "B0079V_post"
# [57] "B0080V_post" "B0081V_post" "B0082V_post" "B0086V_post" "B0088V_post" "B0089V_post" "B0092V_post"
# [64] "B0093V_post" "B0096V_post" "B0098V_post" "B0099V_post" "B0100V_post" "B0102V_post" "B0103V_post"
# [71] "B0104V_post" "B0105V_post" "B0106V_post" "B0107V_post" "B0108V_post" "B0109V_post" "B0110V_post"
# [78] "B0111V_post" "B0112V_post" "B0113V_post" "B0114V_post" "B0115V_post" "B0116V_post" "B0118V_post"
# [85] "B0120V_post" "B0122V_post" "B0123V_post" "B0124V_post" "B0125V_post" "B0126V_post" "B0127V_post"
# [92] "B0129V_post" "B0132V_post" "B0133V_post" "B0135V_post" "B0136V_post" "B0137V_post" "B0138V_post"
# [99] "B0139V_post" "B0140V_post" "B0141V_post" "B0142V_post" "B0143V_post" "B0144V_post" "B0145V_post"
#[106] "B0149V_post" "B0150V_post" "B0151V_post" "B0152V_post" "B0153V_post" "B0154V_post" "B0156V_post"
#[113] "B0157V_post" "B0159V_post" "B0160V_post" "B0162V_post" "B0163V_post" "B0164V_post"
# GENE$ugene_id <- paste(GENE$ugene_id, "X", sep="")
# GENE$ugene_id <- sub('X$', '', GENE$ugene_id)
# View(GENE)
# GENE <- read.csv("ratio.csv", stringsAsFactors=FALSE)
# rownames(LV118pre)=gsub("(B.*V)", "\\1_\\1", rownames(GenoCov120))
# df$ugene_id <- sub('X$', '', df$ugene_id)
LV118pre$ugene_id  <- rownames(LV118pre)
LV118post$ugene_id <- rownames(LV118post)
LV118pre  <- merge(mergedID, LV118pre, by="ugene_id", all.x=TRUE)
LV118post <- merge(mergedID, LV118post, by="ugene_id", all.x=TRUE)
rownames(LV118pre)  <- LV118pre[, 1]
rownames(LV118post) <- LV118post[, 1]
LV118pre  <- LV118pre[, -1]
LV118post <- LV118post[, -1]
rownames(LV118pre)  <- paste(rownames(LV118pre), "X", sep="")
rownames(LV118post) <- paste(rownames(LV118post), "X", sep="")
colnames(LV118pre)  <- sub('_pre$',  '', colnames(LV118pre))
colnames(LV118post) <- sub('_post$', '', colnames(LV118post))
colnames(LV118pre)  <- gsub("(B.*V)", "\\1_\\1", colnames(LV118pre))
colnames(LV118post) <- gsub("(B.*V)", "\\1_\\1", colnames(LV118post))
write.csv(LV118pre,   file="./LV118preX.csv",  row.names=TRUE)
write.csv(LV118post,  file="./LV118postX.csv", row.names=TRUE)
