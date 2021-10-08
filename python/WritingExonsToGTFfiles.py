__author__ = 'layal'

import os
import subprocess
import re
import itertools
import numpy as np
import matplotlib.pyplot as plt
import HTSeq


path = '/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO/GeneAnnotation/Mus_musculus.GRCm38.86.gtf'
root = '/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO/coverage/'

'''
The purpose of this script is to divide up the genes in a GTF file into three files.
The file first_exons.gtf has the 5' exons of each of the genes
last_exons.gtf has the 3' exons
middle_exons has the exons 2..n-1
The script only takes GTF lines about exons and discards other features.
'''

# here we create new empty files and assign them to file object
F = open(root +'first_exons_GRCm38.gtf', 'w')
L = open(root +'last_exons_GRCm38.gtf', 'w')
M = open(root +'middle_exons_GRCm38.gtf', 'w')

'''
Extract a gene id from a GTF line.
The assumption is that the 9th [8] field has a string like this:


gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935"; havana_gene_version "1";
gene_id "ENSMUSG00000102693"; gene_version "1"; transcript_id "ENSMUST00000193812"; transcript_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935"; havana_gene_version "1"; transcript_name "4933401J01Rik-001"; transcript_source "havana"; transcript_biotype "TEC"; havana_transcript "OTTMUST00000127109"; havana_transcript_version "1"; tag "basic"; transcript_support_level "NA";
gene_id "ENSMUSG00000102693"; gene_version "1"; transcript_id "ENSMUST00000193812"; transcript_version "1"; exon_number "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935"; havana_gene_version "1"; transcript_name "4933401J01Rik-001"; transcript_source "havana"; transcript_biotype "TEC"; havana_transcript "OTTMUST00000127109"; havana_transcript_version "1"; exon_id "ENSMUSE00001343744"; exon_version "1"; tag "basic"; transcript_support_level "NA";
gene_id "ENSMUSG00000064842"; gene_version "1"; gene_name "Gm26206"; gene_source "ensembl"; gene_biotype "snRNA";
the first element [0] in the line split by ';' is the id
the second element [1] in the id part split by " " is the id
In this example, we would return the String "ENSMUSG00000102693" , "ENSMUSG00000064842"
'''
def get_geneID(line):
    Ar = line.split('\t')
    id = Ar[8].split(';')[0]
    gene_id = id.split('"')[1]
    return (gene_id)

'''
Plus strand holds all of the GTF lines for genes on the plus strand
 Key is a gene ID, e.g. "ENSMUSG00000064842"
Value is an array of GTF lines
Same as plus_strand but for minus-strand genes
'''
plus_strand = dict()

reverse_strand = dict()

## The following iteration puts the exon definition lines into the plus_strand and minus_strand dictionaries
'''
the lines look like this:
1	havana	gene	3073253	3074322	.	+	.	gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935"; havana_gene_version "1";
1	havana	transcript	3073253	3074322	.	+	.	gene_id "ENSMUSG00000102693"; gene_version "1"; transcript_id "ENSMUST00000193812"; transcript_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935"; havana_gene_version "1"; transcript_name "4933401J01Rik-001"; transcript_source "havana"; transcript_biotype "TEC"; havana_transcript "OTTMUST00000127109"; havana_transcript_version "1"; tag "basic"; transcript_support_level "NA";
'''
with open(path) as gtf:
    for line in gtf:
        #print "line is : ", line , "\n"
        if line.startswith("#!"):
            continue
        else:
            #print "line is (annotation line): ", line , "\n"
            arr = line.split('\t')
            if arr[2] =='exon':
                arr_exon = line.split('\t')
                if arr_exon[6]== '+':
                    gene_idP = get_geneID(line)
                    if gene_idP in plus_strand:
                        plus_strand[gene_idP].append(line)
                        #print "plus strand",plus_strand[gene_idP]
                    else:
                        plus_strand[gene_idP] = [line]
                        #print "plus strand",plus_strand[gene_idP]
                elif arr_exon[6]== '-':
                    gene_idR = get_geneID(line)
                    if gene_idR in reverse_strand:
                        reverse_strand[gene_idR].append(line)
                        #print "reverse", reverse_strand[gene_idR]
                    else:
                        reverse_strand[gene_idR]= [line]
                        #print "reverse", reverse_strand[gene_idR]

## Now write the corresponding data into the three output files separate the exons
for key in plus_strand:
    ar = plus_strand[key]
    ln = len(ar)
    #print "PLUS STRAND ARRAY LINE :" , ar
    #print "first element ::::  ", ar[0]
    #print "length :" , len(ar)
    F.write(ar[0])
    L.write(ar[ln-1])
    for i in range(1, ln-2):
        M.write(ar[i])

for key in reverse_strand:
    ar1 = reverse_strand[key]
    ln = len(ar1)
    F.write(ar1[ln-1])
    L.write(ar1[0])
    for i in range(1, ln-2):
        M.write(ar1[i])