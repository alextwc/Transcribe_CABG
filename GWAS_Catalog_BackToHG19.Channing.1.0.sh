#!/bin/bash
#=============================================================================#
#Title           :GWAS_Catalog_BackToHG19.Channing.1.0.sh                     #
#Description     :The script was written for converting the GWAS_Catalog HG38 #
#                 build back to HG19 build                                    #
#Author          :Dr. Alex Tzuu-Wang Chang                                    #
#Date            :20161018                                                    #
#Version         :1.0 (Channing Cluster Environment Use Only)                 #
#Usage           :qsub -l lx yourScript.sh                                    #
#Require&Notes   :qsub of SGE & BASH shell to use this script                 #
#BASH_Version    :4.1.2(1)-release                                            #
#=============================================================================#
#$ -N convertKGP2RS.Channing.1.0.sh
#$ -M retwc@channing.harvard.edu
#$ -cwd
#$ -l virtual_free=32G
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
echo "## This script takes advantage of \"UCSC LiftOver\" program to reverse the SNP coordinate of HG38 genome    ##"
echo "## assembly back to HG19 genome assembly therefore you have to install \"LiftOver\" onto your cluster first.##"
echo "## To use the LiftOver program you require a UCSC-generated over.chain file as input. Pre-generated over- ##"
echo "## chain files are available from the following web site:                                                 ##"
echo "##        \"wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz\"         ##"
echo "##                                                                                                        ##"
echo "## liftOver binary file can be downloaded from here:                                                      ##"
echo "##              \"wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/lifeOver\"                     ##"
echo "############################################################################################################"
echo "############################################################################################################"
echo
#######################################################################################################
date && echo "Setting the PATH and environmental variables, removing annotated header lines now"
#######################################################################################################
echo
SourceDIR="/path/20161018GWAS_Catalog"
  DestDIR="/path/20161018GWAS_Catalog"
cd $DestDIR

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/lifeOver
ls -alish

echo    
#######################################################################################################
date && echo "Starting downloading the most updated GWAS Catalogue file now"
date && echo "The SNPs with multiple coordinates and empty coordinate will be removed during hg19 coordinate conversion"
date && echo "The output file has been set to BED format to allow directly using in IGViewer"
date && echo "The BED format does not contain the header annotation"
#######################################################################################################
curl -s https://www.ebi.ac.uk/gwas/api/search/downloads/alternative > gwas_catalog_v1.0.1-associations_e86_r2016-10-17.tsv
 awk 'BEGIN {FS="\t"; OFS="\t"} $13 ~ /[0-9]*;;*/ {$13=""} {print}' ./gwas_catalog_v1.0.1-associations_e86_r2016-10-17.tsv\
 | awk '{FS="\t"; OFS="\t"} $13 ~ /[0-9]*xx*/ {$13=""} {print}'\
 | awk '{FS="\t"; OFS="\t"; if(NR>1 && $13!="") print "chr"$12,$13-1,$13,"hg38_chr"$12"_"$13"_"$22,1,".",$8,$11,$15,$21,$22,$24,$26,$27,$28,$31,$35,$37}'\
 | liftOver -bedPlus=4 -tab stdin hg38ToHg19.over.chain.gz stdout unmapped | sort -u > ./gwas_catalog_v1.0.1-associations_e86_r2016-10-17.hg19Conversion.bed

echo    
#######################################################################################################
date && echo "Starting downloading the most updated GWAS Cardiovascular Disease Catalogue file now"
date && echo "The SNPs with multiple coordinates and empty coordinate will be removed during hg19 coordinate conversion"
date && echo "The output file has been set to BED format to allow directly using in IGViewer"
date && echo "The BED format does not contain the header annotation"
#######################################################################################################
curl -s https://www.ebi.ac.uk/gwas/search?query=cardiovascular%20disease/ > ./gwas-association-downloaded_2016-10-18-cardiovascular_disease.tsv
 awk 'BEGIN {FS="\t"; OFS="\t"} $13 ~ /[0-9]*;;*/ {$13=""} {print}' ./gwas-association-downloaded_2016-10-18-cardiovascular_disease.tsv\
 | awk '{FS="\t"; OFS="\t"} $13 ~ /[0-9]*xx*/ {$13=""} {print}'\
 | awk '{FS="\t"; OFS="\t"; if(NR>1 && $13!="") print "chr"$12,$13-1,$13,"hg38_chr"$12"_"$13"_"$22,1,".",$8,$11,$15,$21,$22,$24,$26,$27,$28,$31}'\
 | liftOver -bedPlus=4 -tab stdin hg38ToHg19.over.chain.gz stdout unmapped | sort -u > ./gwas-association-downloaded_2016-10-18-cardiovascular_disease.hg19Conversion.bed

#######################################################################################################
date && echo "Finished the entire AWK script now"
#######################################################################################################
echo
