#!/bin/bash

PathAlign=/home/layal/Solexa/layal/DataOsteo/tophat_out
PathS=/home/layal/PycharmProjects/RSeQC-2.6.2/scripts
PathOut=/home/layal/SVN/DEXseq2Uniprot/results/img


python $PathS/mismatch_profile.py  -i $PathAlign/66604/accepted_hits.bam -l 101 -o $PathOut/DayZeroR1
