__author__ = 'layal'

import os
import subprocess
import re
import socket
import zipfile


def  ( path ):

    #read the zip file in QC which contains fastqc_data.txt
    fp = zipfile.ZipFile(path, 'r')
    inModule= False
    print ("Analysing file ", path)
    for name in fp.namelist(): # iterate through all file names in Zip archive
        if name.find('fastqc_data.txt')>=0:
            for line in fp.read(name).split("\n"): ## split file into lines fp:
             if (line.startswith('#')):
                continue
             if (line.startswith(">>Basic Statistics")):
                inModule=True
                continue
             if (line.startswith(">>END_MODULE")):
                inModule=False
             if (inModule):
                if (line.startswith('Total Sequences')):
                   Arr1 = line.rstrip().split('\t')
                   numberOfReads = Arr1[1]
                if (line.startswith('Sequence length')):
                   Arr = line.rstrip().split('\t')
                   length = Arr[1]
    return (numberOfReads, length)