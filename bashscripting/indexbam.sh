#!/bin/bash

# The purpose of this script is to index the bam aligned files
# Input: the aligned files by tophat
# samtools version 1.3.1
# Output: index .bam.bai files

#path to where the data lives
BASE=/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO

# the input bam files
tophat_out=$BASE/sorted_Bam_GRCm38


if [ ! -d ${tophat_out} ]; then
    mkdir ${tophat_out}
fi

#index the BAM files for each sample according to the position of reference genome :


for F in $tophat_out/*.bam; 
do
    i=$(basename $F); ## is is the name of the file, e.g., 666004_sorted.bam
    
    echo "Calling samtool for INDEXING on $i "
        
   samtools index  -b ${tophat_out}/${i} ${tophat_out}/${i}.bai
done



