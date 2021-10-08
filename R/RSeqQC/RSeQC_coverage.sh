#!/bin/bash

PathAlign=/home/layal/Solexa/layal/DataOsteo/tophat_out
PathRef=/home/layal/Solexa/layal/DataOsteo
Path=/home/layal/PycharmProjects/RSeQC-2.6.2/scripts

#geneBody_coverage.py -r $Path/mm10.HouseKeepingGene.bed -i $Path/bam_path.txt -l 50 -o output


python $Path/geneBody_coverage.py -r $PathRef/mm10.HouseKeepingGene.bed -i $PathAlign/66604/accepted_hits.bam,$PathAlign/66605/accepted_hits.bam,$PathAlign/66606/accepted_hits.bam,$PathAlign/66607/accepted_hits.bam,$PathAlign/66608/accepted_hits.bam,$PathAlign/66609/accepted_hits.bam,$PathAlign/66610/accepted_hits.bam,$PathAlign/66611/accepted_hits.bam,$PathAlign/66613/accepted_hits.bam,$PathAlign/66614/accepted_hits.bam,$PathAlign/66615/accepted_hits.bam,$PathAlign/66616/accepted_hits.bam -l 101 -o /home/layal/SVN/DEXseq2Uniprot/results/output


