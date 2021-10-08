#!/bin/bash

#the script calls out Trimmomatic to perform trimming on the RNAseq reads

#path to where the data lives in solexa
BASE='/home/layal/OSTEO'

TRIMMOMATIC='/home/layal/Trimmomatic-0.36'

## change path for peter's computer
if [ $HOSTNAME = 'bioinf-hpo' ]; then
    BASE=/home/peter/data/rnaseq
    TRIMMOMATIC='/home/peter/bin/Trimmomatic-0.32'
fi

pathRNA=${BASE}/RNAseq_raw

pathTrimRNAout=${BASE}/RNA_afterTrimming

## make folder of the RNAseq the output of the trimming
if [ ! -d ${pathTrimRNAout} ];
then
    mkdir ${pathTrimRNAout}
fi

for Dir in $(find $pathRNA -mindepth 1 -maxdepth 1 -type d ); 
do
    i=$(basename $Dir); ## is is the name of the experiment, e.g., 666004
    echo "Calling Trimmomatic on $Dir for experiment $i"
    mkdir $pathTrimRNAout/${i}
    java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads 10 -phred33 -trimlog  $pathTrimRNAout/${i}.log   $pathRNA/${i}/${i}_1.fq.gz  $pathRNA/${i}/${i}_2.fq.gz  $pathTrimRNAout/${i}/${i}_1_paired_output.fq.gz  $pathTrimRNAout/${i}/${i}_1_unpaired_output.fq.gz  $pathTrimRNAout/${i}/${i}_2_paired_output.fq.gz  $pathTrimRNAout/${i}/${i}_2_unpaired_output.fq.gz  ILLUMINACLIP:${TRIMMOMATIC}/adapters/TruSeq3-PE.fa:2:30:10:4:true LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36
done


# ========================== For easier further process we want to move unpaired RNA seq and rename "*_paired_output.fq " as "*.fq" ================================
# make directory for the unpaired RNA (the output from Trimming) where we want to move them. 
pathunpairedRNA=${BASE}/RNA_afterTrimming_unpaired

if [ ! -d ${pathunpairedRNA} ];
then
    mkdir ${pathunpairedRNA}
fi



for Dir in $(find $pathTrimRNAout -mindepth 1 -maxdepth 1 -type d )
do
    i=$(basename $Dir); ## is is the name of the experiment, e.g., 666004
    echo "movinig unpaired RNA from  $Dir for experiment $i to $pathunpairedRNA "
    mkdir $pathunpairedRNA/${i}
    find $pathTrimRNAout/${i} -type f -name "*_unpaired_output.fq.gz" -exec mv {} $pathunpairedRNA/${i} \;
    echo "rename (_paired_output.fq.gz) files for easier uses as this (*_1.fq.gz) " 
    find $pathTrimRNAout/${i} -type f -name "*1_paired_output.fq.gz" -exec mv {} $pathTrimRNAout/${i}/"${i}_1.fq.gz" \;
    find $pathTrimRNAout/${i} -type f -name "*2_paired_output.fq.gz" -exec mv {} $pathTrimRNAout/${i}/"${i}_2.fq.gz" \;
done
    





