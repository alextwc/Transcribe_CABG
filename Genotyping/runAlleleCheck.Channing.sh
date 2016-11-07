#!/bin/bash
#$ -N runAlleleCheck.Channing.sh
#$ -M retwc@channing.harvard.edu
#$ -cwd
#$ -l virtual_free=32G
#$ -o log.$JOB_NAME.$JOB_ID.txt 
#$ -e err.$JOB_NAME.$JOB_ID.txt
#$ -S /bin/bash

: <<ScriptBackground
This BASH script is for the automation of 
"running John Ziniti's allele_check.sh".
Script Author: Dr. Alex Tzuu-Wang Chang 
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
# For i  in  $(ls  $SourceDIR)
# do 
#         if [ -d $i ]; then
#              break
#         fi 
#         NewFilename=$( ls $i | awk -F-- 'BEGIN{OFS="-"}{print $3, $4, $2, $1, $5}')
#         $i=$NewFilename 
# done
# cd "/path/analyses/retwc/Temp/temp"
# tar -cvfz 35subjects.tar.gz /path/analyses/retwc/Temp/temp/20140623_35Subjects
# TMP_file="/tmp/tmp.$USER.$$"

dtime=$(date +'%F_%H%M')
echo "The current job name is $JOB_NAME"
echo "The current time is $dtime"
echo "The current shell ID is $$"
echo "The JOB_ID of this submitted shell script (BASH script) is $JOB_ID"

cd /path/analyses/retwc/Omni2-5_round_2
( source ./allele_check.sh tmp/munged.bim ) >& runAlleleCheckLog.txt 
