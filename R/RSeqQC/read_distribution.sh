#!/bin/bash

PathAlign=/home/layal/Solexa/layal/DataOsteo/tophat_out
Path=/home/layal/Solexa/layal/DataOsteo
PathS=/home/layal/PycharmProjects/RSeQC-2.6.2/scripts

python $PathS/read_distribution.py -i $PathAlign/66604/accepted_hits.bam -r $Path/mm10_RefSeq.bed
