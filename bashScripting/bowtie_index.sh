#!/bin/bash

# build the index  of GRCm38, the fasta file must be unzipped
# option -f  The reference input files are FASTA files (usually having extension .fa, .mfa, .fna or similar).
#Run in solexa : /home//layal/OSTEO/Bowtie2_index
# after we get the FASTA file from Ensemble Mus_musculus.GRCm38.dna.toplevel.fa.gz (I think we need to call it as BASENAME.fa )
# we unzip it by : gunzip -k Mus_musculus.GRCm38.dna.toplevel.fa.gz.
# then we run the following command

bowtie2-build -f /home/layal/OSTEO/Bowtie2_index/Mus_musculus.GRCm38.fa /home/layal/OSTEO/Bowtie2_index/Mus_musculus.GRCm38
