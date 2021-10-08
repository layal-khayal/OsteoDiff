#!/bin/bash

PathAlign=/home/layal/Solexa/layal/DataOsteo/tophat_out
Path=/home/layal/PycharmProjects/RSeQC-2.6.2/scripts


python $Path/read_duplication.py -i $PathAlign/66604/accepted_hits.bam -o /home/layal/SVN/DEXseq2Uniprot/results/DayZero
