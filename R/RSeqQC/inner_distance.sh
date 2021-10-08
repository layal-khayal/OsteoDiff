#!/bin/bash

PathAlign=/home/layal/Solexa/layal/DataOsteo/tophat_out
PathRef=/home/layal/Solexa/layal/DataOsteo
Path=/home/layal/PycharmProjects/RSeQC-2.6.2/scripts
PathOut=/home/layal/SVN/DEXseq2Uniprot/results

python $Path/inner_distance.py -i $PathAlign/66604/accepted_hits.bam -o $PathOut/66604_output -r $PathRef/mm10_RefSeq.bed

