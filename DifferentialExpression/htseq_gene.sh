#!/bin/bash
source /proj/sadevs/sadev01/python/2.7.3/bin/activate
python -m HTSeq.scripts.count -f bam -r name -s reverse $1 $2 > $3