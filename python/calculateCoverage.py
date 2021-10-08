__author__ = 'layal'

import os
import subprocess
import re
import itertools
import numpy as np
import matplotlib.pyplot as plt
import HTSeq
import collections

'''
Purpose of the script:
To calculate the coverage of the exons for some BAM file.
'''
# This function calculate the coverage of one feature with is exons from
# an alignment Bam file
#  and the new gtf file we generate for the first exons , last one and the middle ones

def calculateCoverageOfOneBamFile(gtf_file,Alignment, filename_cv):
    C = open(filename_cv, 'w')

    # we create a GennomicArray data structure to store and retrieve information associated with a genomic position or genomic interval.
    # we instruct the GenomicArray to add chromosome vectors as needed, by specifiyng "auto".
    # the stranded is false it is not relevant the direction
    # typecode='i' because we are saving intager in the genomic Array .
    # This will store intager of Genomic intervals

    print "Asign Genomic array"
    cvg= HTSeq.GenomicArray("auto" , stranded=False, typecode = 'i')

    #We are counting the coverage by adding one to the internval of the read which is aligned, cvg is a GenomicArray object

    print "looping over the aligned object to set values in the Genomic array"
    for align in Alignment:
        if align.aligned:
            cvg[align.iv] +=1


    # now we have the interval of the aligned read assigned by adding one each time,
    # gtf file is an class object from GFF_Reader
    MeanCoverage = []
    print "looping through the features in each gtf file"
    for feature in gtf_file:
        if feature.type == "exon":
            mv = np.mean(list(cvg[ feature.iv ]))
            C.write(feature.name + "\t" + format(str(mv) + "\n"))
            MeanCoverage.append(mv)

    return(MeanCoverage)






