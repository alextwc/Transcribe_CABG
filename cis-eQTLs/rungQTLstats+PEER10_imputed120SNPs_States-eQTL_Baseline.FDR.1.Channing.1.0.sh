#!/bin/bash
#=============================================================================#
#Title           :rungQTLstats+PEER10_imputed120SNPs_States-eQTL_Baseline.FDR #
#                 .1.Channing.1.0.sh                                          #
#Description     :This script will run a R script                             #
#Author          :Dr. Alex Tzuu-Wang Chang                                    #
#Date            :2016-09-27                                                  #
#Version         :1.0                                                         #
#Usage           :qsub -l lx ./yourScript.Channing.x.x.sh                     #
#Require&Notes   :qsub of SGE & BASH shell to use this script                 #
#BASH_Version    :4.1.2(1)-release                                            #
#=============================================================================#
#$ -N rungQTLstats+PEER10_imputed120SNPs_States-eQTL_Baseline.FDR.1.Channing.1.0.sh
#$ -M retwc@channing.harvard.edu
#$ -cwd
#$ -l virtual_free=92G
#$ -l m_core=9
#$ -o log.$JOB_NAME.$JOB_ID.txt 
#$ -e err.$JOB_NAME.$JOB_ID.txt
#$ -S /bin/bash
: <<ScriptBackground
This BASH script is for the automation of 
"Running R script in background with unattended mode".
ScriptBackground

# . /proj/sadevs/sadev01/settings/pipeline_env.sh
# . /proj/sadevs/sadev01/python/2.7.3/bin/activate
# export PATH=/local/bin/:$PATH
# export TMPDIR=/udd/retwc/Temp
# export R_BZIPCMD="/usr/bin/bzip2"
# export R_INCLUDE_DIR="/app/R-3.0.2@i86-rhel6.0/lib64/R/include"
# export R_PAPERSIZE="a4"
# export R_PRINTCMD="lpr"
# export R_SHARE_DIR="/app/R-3.0.2@i86-rhel6.0/lib64/R/share"
# export R_UNZIPCMD="/usr/bin/unzip"
# export R_ZIPCMD="/usr/bin/zip"
# export R_SYSTEM_ABI="linux,gcc,gxx,gfortran,?"
# export R_RD4PDF="times,hyper"
# export R_PDFVIEWER="/usr/bin/xdg-open"
# export R_GZIPCMD="/bin/gzip"
# export R_DOC_DIR="/app/R-3.0.2@i86-rhel6.0/lib64/R/doc"
# export R_HOME="/app/R-3.0.2@i86-rhel6.0/lib64/R"
# export R_PLATFORM="x86_64-unknown-linux-gnu"
# export R_TEXI2DVICMD="texi2dvi"
# echo 'R_LIBS_USER="~/R/library"' >  $HOME/.Renviron
# export R_LIBS_USER="/home/tac15/R/library"
# For i  in  $(ls  $SourceDIR)
# do 
#         if [ -d $i ]; then
#              break
#         fi 
#         NewFilename=$( ls $i | awk -F-- 'BEGIN{OFS="-"}{print $3, $4, $2, $1, $5}')
#         $i=$NewFilename 
# done
# cd "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/Temp/temp"
# tar -cvfz 35subjects.tar.gz /proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/Temp/temp/20140623_35Subjects
# TMP_file="/tmp/tmp.$USER.$$"
# dtime=$(date +'%F_%H%M')
# SCRIPT_SHELL=${SHELL}
# DATE=$(date +'%Y-%m-%d')
# YEAR=$(date +'%Y')

dtime=`date +'%F_%H%M'`
echo "The current job name is $JOB_NAME"
echo "The current time is $dtime"
echo "The current shell ID is $$"
echo "The JOB_ID of this submitted shell script (BASH script) is $JOB_ID"
echo

date && echo "Setting paths for Source and Destination"
SourceDIR="/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160923reRUN01"
  DestDIR="/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20160923reRUN01"
echo

date && echo "Switching to Working Directory then creating list"
cd $DestDIR
: > ./$JOB_NAME.txt
echo

date && echo "Starting running R script"
( /app/R-3.3.1@i86-rhel6.0/bin/Rscript ./gQTLstats+PEER10_imputed120SNPs_States-eQTL_Baseline.FDR.1.Channing.1.0.R ./InputFiles/LV135preX.csv ./InputFiles/LV133postX.csv ./InputFiles/CovariancesAll.csv) >& run.$dtime.$$.$JOB_NAME.log
echo

date && echo "Writing the list of files"
echo >> $JOB_NAME.txt
ls -alish >> $JOB_NAME.txt
echo

date && echo "Finished"
