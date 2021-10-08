#!/bin/bash

## This script processes the raw data with FastQC v0.11.5
## It puts the output into a new directory called FASTQC_OUTPUT

#path to where the data lives in solexa
BASE=/home/layal
FASTQC=${BASE}/FastQC/fastqc


## change path for peter's computer
if [ $HOSTNAME = 'bioinf-hpo' ]; then
    BASE=/home/peter/data/rnaseq
    FASTQC=/home/peter/bin/FastQC/fastqc
fi

# Path to the raw FASTQ files
pathRNA=${BASE}/OSTEO/RNAseq_raw

# Path to where we write the output data
out_path=${BASE}/OSTEO/FASTQC_OUTPUT

## Create the directory. We will write FASTQC output to a subdirectory called QC
if [ ! -d $out_path ]; then
    echo "Making directory $out_path"
    mkdir $out_path
fi

if [ ! -d ${out_path}/QC ]; then
    echo "Making directory ${out_path}/QC"
    mkdir ${out_path}/QC
fi


for Dir in $(find $pathRNA -mindepth 1 -maxdepth 1 -type d ); 
do
    i=$(basename $Dir); ## is is the name of the experiment, e.g., 666004
    echo "Calling FASTQC on $Dir for experiment $i"
    ## Call FASTQC for forward and reverse strands
    $FASTQC ${Dir}/${i}_1.fq.gz -o $out_path/QC/
    $FASTQC ${Dir}/${i}_2.fq.gz -o $out_path/QC/
done





