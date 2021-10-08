#!/bin/bash

# The purpose of this script is to run tophat2 on our Osteoblast data (align the reads)
# Input: the trimmed RNAseq files 666004, 666005, ..., 66616 (note: 66612 was skipped)
# Input: bowtie2 index files for GRCm38
# Bowtie2 version 2.2.4 and Topghat is 2.0.13
# gene annotation file version GRCm38 release 86
# Output: Mapped reads at $Align_Tophat_GRCm38


#path to where the data lives $on Solexa
BASE=/home/layal/OSTEO

#notice the gtf file must be unzip
geneAnnot=$BASE/GeneAnnotation/Mus_musculus.GRCm38.86.gtf
genomeIndex=$BASE/Bowtie2_index/Mus_musculus.GRCm38
tophat_out=$BASE/Align_Tophat_GRCm38
RNAseq_data=$BASE/RNA_afterTrimming
TOPHAT=/home/layal/tophat/tophat2



#must be changed to annotation Mus_musculus.GRCm38.86.gtf and Mus_musculus.GRCm38 Bowtie index
## change path for peter's computer
if [ $HOSTNAME = 'bioinf-hpo' ]; then
    BASE='/home/peter/data/rnaseq'
    geneAnnot="${BASE}/mm10/genes.gtf"
    genomeIndex="${BASE}/mm10/Bowtie2Index/genome"
    RNAseq_data="${BASE}/RNA_afterTrimming"
    TOPHAT='/home/peter/bin/tophat-2.1.1.Linux_x86_64/tophat2'
fi

if [ ! -d ${tophat_out} ]; then
    mkdir ${tophat_out}
fi




#Map the reads for each sample to the reference genome:


for Dir in $(find $RNAseq_data -mindepth 1 -maxdepth 1 -type d ); 
do
    i=$(basename $Dir); ## is is the name of the experiment, e.g., 666004
    echo "Calling tophat on $Dir for experiment $i"
    mkdir $tophat_out/${i}
    ${TOPHAT} -p 30 -G $geneAnnot --b2-sensitive --keep-fasta-order -o $tophat_out/${i} $genomeIndex $RNAseq_data/${i}/${i}_1.fq.gz $RNAseq_data/${i}/${i}_2.fq.gz
done
