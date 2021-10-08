__author__ = 'layal'

import os
import subprocess
import re
import sys
import socket
import zipfile



def BasicStatistic ( path ):

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

def AnalyzeDataset (path , outputname):
    '''Analyse a set of Fastqc output files and extract the mean
    quality values
    '''
    print ("Path" , path)
    F = open(outputname, 'w')
    F.write(' \n The length and number of read before trimming\n\n')
    # os.walk function return dirpath is a string, the path to the directory.
    #  dirnames is a list of the names of the subdirectories in dirpath . filenames is a list of the names of the non-directory files in dirpath.
    for dirPath, dirNames, fileList in os.walk(path):
        #print("dirPath: " , dirPath)
        #print("dirNames: ", dirNames , " fileList:" , fileList)
        for index, item in enumerate(fileList):
            #print("item :", item)
            if item.endswith("zip"):
                print("item :", item)
                dpath= dirPath + '/' + item
                print ("PATH ::::", dpath)
                pos = dpath.rfind('/')
                sample = dpath[pos+1:pos+6]
                print ("Analyzing directory: " , sample)
                numReads, leng = BasicStatistic(dpath)
                print ('number of reads = ' , numReads ,'\n reads length = ' , leng )
                F.write('Sample : '+ str(sample) + '\t\t Reads number : '+ str(numReads) +'\t\t Reads length :' + str(leng) + '\n' )
            enumerate


## Path to the directory with the FastQC files. before trimming

results_QC_before = '/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/FASTQC_OUTPUT/QC'
AnalyzeDataset(results_QC_before, "basicStatistics_before_trimming.txt")



results_QC_after = '/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/FASTQC_OUTPUT/QC_aftertrimming'

AnalyzeDataset(results_QC_after, "basicStatistics_after_trimming.txt")
