#!/bin/bash

# The purpose of this script is to sort the aligned files according to position
# Input: the aligned files by tophat
# samtools version 1.3.1
# Output: soted Bamfiles

#path to where the data lives
BASE=/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO

# the input bam files
tophat_out=$BASE/Align_Tophat_GRCm38
Sorted_Bam=$BASE/sorted_Bam_GRCm38

if [ ! -d ${Sorted_Bam} ]; then
    mkdir ${Sorted_Bam}
fi

#Sort the BAM files for each sample according to the position of reference genome :


for Dir in $(find $tophat_out -mindepth 1 -maxdepth 1 -type d ); 
do
    i=$(basename $Dir); ## is is the name of the experiment, e.g., 666004
    echo "Calling samtool for SORTING on $Dir for experiment $i"
        
   samtools sort  -o ${Sorted_Bam}/${i}_sorted.bam ${tophat_out}/${i}/accepted_hits.bam
done




