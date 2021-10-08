__author__ = 'layal'

import os
import subprocess
import re
import socket
import zipfile



'''
The following function searches for the part of the
fastqc_data.txt file (which contains quality scores)
and uses them to calculate the mean quality across all
bases in the BAM file.
>>Per sequence quality scores	pass
#Quality	Count
2	1290108.0
3	314493.0
4	170961.0
(5-39)
40	160.0
>>END_MODULE

Uncomment final command to analyse data after trimming

'''

## Path to FastQC output files (before trimming)
#results_QC_before = '/home/layal/Solexa/layal/OSTEO/FASTQC_OUTPUT/QC'
results_QC_before ='/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/FASTQC_OUTPUT/QC'

## Change path on Peter's computer
hname = socket.gethostname()
if hname == 'bioinf-hpo':
    results_QC_before = '/home/peter/data/rnaseq/FASTQC_OUTPUT/QC'



def QualityScore( path ):
    '''Calculate the mean quality score of the sequence reads in the
    FASTQ file at path. Return the mean quality
    Path is a Zip archive made by FastQC
    '''
    archive = zipfile.ZipFile(path, 'r')
    #  one of the files is 'fastqc_data.txt'
    total_base=0
    sumqual=0
    meanqual=0
    inModule= False
    for name in archive.namelist(): # iterate through all file names in Zip archive
        if name.find('fastqc_data.txt')>=0:
            for line in archive.read(name).split("\n"): ## split file into lines
                if (line.startswith('#')):
                    continue
                if (line.startswith(">>Per sequence quality scores")):
                    inModule=True
                    continue
                if (line.startswith(">>END_MODULE")):
                    inModule=False
                    continue
                if (inModule):
                    ar = line.rstrip().split('\t')
                    phred = int(ar[0])
                    num = float(ar[1]) ## Number of observations, note we cast to int here
                    #print phred,num
                    sumqual+=phred*num
                    total_base+=num
                    meanqual=sumqual/total_base
            return (meanqual)
             

def AnalyzeDataset (path , outputname, Av_mean_output):
    '''Analyse a set of Fastqc output files and extract the mean
    quality values
    '''

    meanqual_samples = {}
    print (path)
    F = open(outputname, 'w')
    F.write(outputname[0:-3] +'\n\n')
    H = open(Av_mean_output, 'w')
    H.write(Av_mean_output[0:-3]+'\n\n')

    # os.walk function return dirpath is a string, the path to the directory.
    #  dirnames is a list of the names of the subdirectories in dirpath . filenames is a list of the names of the non-directory files in dirpath.
    for dirPath, dirNames, fileList in os.walk(path):
        #print("dirPath: " , dirPath)
        #print("dirNames: ", dirNames , " fileList:" , fileList)
        for index, item in enumerate(fileList):
            if item.endswith("zip"):   # pick out the Zip archives from FastQC
                dpath= dirPath + '/' + item
                print ("PATH ::::", dpath)
                pos = dpath.rfind('/')
                sample = dpath[pos+1:pos+6]
                print ("Analyzing directory: " , sample)
                mqual = QualityScore(dpath)
                print ('mean quality = ' , mqual)
                if sample in meanqual_samples:
                    meanqual_samples[sample] += mqual
                else: meanqual_samples[sample] = mqual
                # file contains all the mean quality of all samples
                F.write('Sample: '+ str(item) + '\t\t\t mean quality: '+ str(round(mqual,3))+'\n' )
            enumerate
    # to calculate the average mean quality forward and reverse RNAseq


    for key, value in meanqual_samples.iteritems():
        meanqual_samples[key] = value / 2
        print (key, meanqual_samples[key])
        H.write('Sample: '+ str(key) + '\t\t\t mean quality: '+ str(round(meanqual_samples[key],3))+'\n' )


AnalyzeDataset(results_QC_before, "meanquality_before_trimming.txt", "averrage_meanqual_before_trimming.txt")

## Path to the directory with the FastQC files after trimming.

results_QC_after = '/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/FASTQC_OUTPUT/QC_aftertrimming'

AnalyzeDataset(results_QC_after, "meanquality_after_trimming.txt", "averrage_meanqual_after_trimming.txt")


