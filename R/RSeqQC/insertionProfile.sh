#!/bin/bash

PathAlign=/home/layal/Solexa/layal/DataOsteo/tophat_out
Path=/home/layal/PycharmProjects/RSeQC-2.6.2/scripts


python $Path/insertion_profile.py -s "PE" -i $PathAlign/66604/accepted_hits.bam -o /home/layal/SVN/DEXseq2Uniprot/results/DayZeroR1
