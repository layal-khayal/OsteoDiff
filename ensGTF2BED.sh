#!/bin/bash

# The purpose of this script is to convert the GTF from Ensemble to BED file

#Annot=/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO/GeneAnnotation
Annot=/home/layal/OSTEO/GeneAnnotation

awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $Annot/Mus_musculus.GRCm38.86.gtf | gtf2bed - > $Annot/Mus_musculus.GRCm38.86.bed
