#!/bin/bash
source /proj/sadevs/sadev01/python/2.7.3/bin/activate
python /udd/rerpc/R/x86_64-unknown-linux-gnu-library/3.2/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos -s reverse $2 $1 $3


