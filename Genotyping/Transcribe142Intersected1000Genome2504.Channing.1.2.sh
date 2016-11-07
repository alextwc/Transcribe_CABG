#!/bin/bash
#=============================================================================#
#Title           :Transcribe142Intersected1000Genome2504.Channing.1.2.sh      #
#Description     :This script was written for intersecting the SNPs between   #
#                 1000Genome Phase3 dataset and our Transcribe142 genotype    #
#                 dataset and replacing our Transcribe142 SNP-ID format with  #
#                 SNP-ID format of 1000 Genome Phase3 dataset.                #
#Author          :Dr. Alex Tzuu-Wang Chang                                    #
#Date            :20160829                                                    #
#Version         :1.2 (Fixed the bug of output file format to tab-delimited)  #
#Usage           :qsub -l lx yourScript.sh                                    #
#Require&Notes   :qsub of SGE & BASH shell to use this script                 #
#BASH_Version    :4.1.2(1)-release                                            #
#=============================================================================#
#$ -N Transcribe142Intersected1000Genome2504.Channing.1.2.sh
#$ -M retwc@channing.harvard.edu
#$ -cwd
#$ -l virtual_free=92G
#$ -o log.$JOB_NAME.$JOB_ID.txt 
#$ -e err.$JOB_NAME.$JOB_ID.txt
#$ -S /bin/bash

: <<ScriptBackground
This BASH script is for the automation of 
"Running BASH script in background with unattended mode".
ScriptBackground

dtime=`date +'%F_%H%M'`
echo "The current job name is $JOB_NAME"
echo "The current time is $dtime"
echo "The current shell ID is $$"
echo "The JOB_ID of this submitted shell script (BASH script) is $JOB_ID"
echo

#######################################################################################################
date && echo "Starting the AWK script now"
#######################################################################################################
echo
echo "############################################################################################################"
echo "############################################################################################################"
echo "## You want to change the SNP_ID of your own genotyping dataset to the SNP_ID format of 1000 GenomePhase3 ##"
echo "## genotyping dataset when both SNPs are at the same SNP locus                                            ##"
echo "## The input files require the following conditions:                                                      ##"
echo "## (1). Both input files have to be VCF file format                                                       ##"
echo "## (2). The additional annotated head lines have to be removed. Only allow the folloiwng header line:     ##"
echo "## #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SID1  SID2  SID3  SID4  SID5  SID6  SID7  ..... ##"
echo "##                                (The SID1 .....SIDN means the subject-ID)                               ##"
echo "## (3). The algorithm of this script is to take the SNPs intersection first based on their SNP loci then  ##"
echo "## start changing the SNP-ID to comply with the RS naming system of 1000 Genome genotyping dataset format ##"
echo "## (4). You put 1000 genome reference file as the first input file and the file you want to re-write (    ##"
echo "## which is transcribe142_chr1.vcf) as the second input file.                                             ##"
echo "############################################################################################################"
echo "############################################################################################################"
echo
#######################################################################################################
date && echo "Setting the PATH and environmental variables, removing annotated header lines now"
#######################################################################################################
echo
SourceDIR="/path/analyses/retwc/20160829PCA"
  DestDIR="/path/analyses/retwc/20160829PCA"
cd $DestDIR
# gunzip ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
( sed '1,252d' ./ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf ) > ALL.chr1.phase3.20130502.temp1.vcf
mv -f ALL.chr1.phase3.20130502.temp1.vcf ALL.chr1.phase3.20130502.vcf
( sed '1,6d' ./Transcribe142_Chr01_nonimputed.vcf ) > Transcribe142_Chr01_nonimputed.temp1.vcf
# mv -f Transcribe142_Chr01_nonimputed.temp1.vcf Transcribe142_Chr01_nonimputed.vcf
NR=`cat ./ALL.chr1.phase3.20130502.vcf | wc -l`
NC=`sed -n '1p' ./ALL.chr1.phase3.20130502.vcf | wc -w`
date && echo "The 1000 Genome reference file contains $NR rows"
date && echo "The 1000 Genome reference file contains $NC columns"
NR=`cat ./Transcribe142_Chr01_nonimputed.temp1.vcf | wc -l`
NC=`sed -n '1p' ./Transcribe142_Chr01_nonimputed.temp1.vcf | wc -w`
date && echo "Your own genotyping dataset contains $NR rows"
date && echo "Your own genotyping dataset contains $NC columns"
date && echo "Finished removing the additional annotated header lines now"
echo
#######################################################################################################
date && echo "Intersecting the common shared SNPs by their SNPs loci now"
#######################################################################################################
# ( awk 'NR==FNR {a[$2] = 0; next} {if($2 in a){print}}' ./Transcribe142_Chr01_nonimputed.vcf ./ALL.chr1.phase3.20130502.vcf ) > test1.vcf
  ( awk 'NR==FNR {n[$2] = $0;  next} {if($2 in n){print}}' ./Transcribe142_Chr01_nonimputed.temp1.vcf ./ALL.chr1.phase3.20130502.vcf ) > commShared_chr1.vcf
# The SNP_ID format of final output belong to 1000 Genome format becuase we put ALL.chr1.phase3.20130502.vcf as the second input argument
NR=`cat ./commShared_chr1.vcf | wc -l`
NC=`sed -n '1p' ./commShared_chr1.vcf | wc -w`
date && echo "The intersected SNPs between two files contains $NR rows"
date && echo "The intersected SNPs between two files contains $NC columns"
date && echo "Finished extracting the common shared SNPs between two files now"
echo    
#######################################################################################################
date && echo "Changing the SNP-ID of our smaller genotyping dataset to RS-naming system of 1000 Genome dataset now"
#######################################################################################################
(awk 'BEGIN {OFS="\t"} NR==FNR {a[FNR]=$2;b[FNR]=$3;next} { for (i in b){if($2==a[i]){$3=b[i];print $0}}}' ./commShared_chr1.vcf ./Transcribe142_Chr01_nonimputed.temp1.vcf) > rsSNPs_Transcribe142_Chr01_nonimputed.vcf
NR=`cat ./rsSNPs_Transcribe142_Chr01_nonimputed.vcf | wc -l`
NC=`sed -n '1p' ./rsSNPs_Transcribe142_Chr01_nonimputed.vcf | wc -w`
date && echo "The final intersected SNPs between two files contains $NR rows"
date && echo "The final intersected SNPs between two files contains $NC columns"
date && echo "Finished changing the SNP-ID format with 1000 Genome RS-naming system now"
echo   
#######################################################################################################
date && echo "Finished the entire AWK script now"
#######################################################################################################
echo
