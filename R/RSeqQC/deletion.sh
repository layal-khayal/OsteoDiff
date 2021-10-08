#!/bin/bash

PathAlign=/home/layal/Solexa/layal/DataOsteo/tophat_out
Path=/home/layal/PycharmProjects/RSeQC-2.6.2/scripts

#deletion_profile.py -i $Path/mm10.HouseKeepingGene.bed -i $Path/bam_path.txt -l 50 -o output

python $Path/deletion_profile.py  -i $PathAlign/66604/accepted_hits.bam -l 101 -o /home/layal/SVN/DEXseq2Uniprot/results/DayZeroR1
