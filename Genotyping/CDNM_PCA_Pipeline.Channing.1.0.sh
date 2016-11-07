#!/bin/bash
#=============================================================================#
#Title           :CDNM_PCA_Pipeline.Channing.1.0.sh                           #
#Description     :Run CDNM_PCA_Trans81                                        #
#Author          :Dr. Alex Tzuu-Wang Chang                                    #
#Date            :20161004                                                    #
#Version         :1.0 (Channing Cluster Environment Use Only)                 #
#Usage           :qsub -l lx yourScript.sh                                    #
#Require&Notes   :qsub of SGE & BASH shell to use this script                 #
#BASH_Version    :4.1.2(1)-release                                            #
#=============================================================================#
#$ -N CDNM_PCA_Pipeline.Channing.1.0.sh
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

qrsh -l lx6 -l mem_free=92G
cd path/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004

plink2 --noweb \
       --file ./Transcribe142_AllChr22_imputed_MAF0.05passed \
       --remove ./RemoveSubjects.txt \
       --out  ./Transcribe081_AllChr22_imputed_MAF0.05passed \
       --recode

plink2 --noweb \
       --make-bed \
       --out  ./Transcribe081_AllChr22_imputed_MAF0.15passed \
       --file ./Transcribe081_AllChr22_imputed_MAF0.05passed \
       --mind 1.0 \
       --geno 1.0 \
       --maf 0.15 \
       --hardy \
       --hwe 0.000001
       
plink2 --noweb \
       --bfile ./Transcribe081_AllChr22_imputed_MAF0.15passed \
       --out   ./Transcribe081_AllChr22_imputed_MAF0.15passed \
       --recode
              
# Becuase the PLINK2 can't create folders therefore we have to run the following commands to build those sub-folders first:
mkdir path/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source
   cd path/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source
for chr in `seq 25`; do mkdir chr${chr}; done
cd ..
mkdir PLINK
   cd PLINK
mv -f ../Transcribe081_AllChr22_imputed_MAF0.15passed.* .
# Now we can start split the All22 chromosomes into each chromosome
date && echo "Splitting BED/BIM/FAM files into 22 chromosome ${chr} working directory now"
for chr in `seq 22`; do plink --nonfounders --allow-no-sex --noweb --bfile Transcribe081_AllChr22_imputed_MAF0.15passed --chr ${chr} --make-bed --out /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr${chr}/PCA_chr${chr}; done
# copy chr22*.* to chr23, chr24 and chr25
for chr in {23,24,25}
do
    cd /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr${chr}
    cp /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr22/PCA_chr22.bed ./PCA_chr${chr}.bed
    cp /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr22/PCA_chr22.bim ./PCA_chr${chr}.bim
    cp /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr22/PCA_chr22.fam ./PCA_chr${chr}.fam
    cp /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr22/PCA_chr22.log ./PCA_chr${chr}.log
done
# Make PED/MAP files from BED/BIM/FAM files
echo
date && echo "Converting MAF=0.15 filtered BED/BIM/FAM files of chromosome ${chr} to regenerate filtered PED/MAP files now"
for chr in `seq 25`; do plink --noweb --bfile /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr${chr}/PCA_chr${chr} --out /udd/retwc/TRANSCRIBE/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/chr${chr}/PCA_chr${chr} --recode; done
cd ..
mkdir path/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Results

#run PCA#
# qrsh -l lx6 -l mem_free=92G
source /proj/relibs/relib00/python/2.7.3_x86_64_CentOS6.5/bin/activate
cd path/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Results
cdnm_pca_pipeline --parallel=100 path/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/
cdnm_pca_pipeline --out-base=Transcribe081_AllChr22_imputed_MAF0.15passed path/data/rnaseq/TRANSCRiBE142Genotypes/imputed/20161004/PCA_Source/
make
# Finished running PCA#
