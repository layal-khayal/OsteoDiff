__author__ = 'layal'
import HTSeq
import itertools
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import socket



pathout= "/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO/RNA_afterTrimming/"



def reads_qual (fq):
    readlen = 101
    qualsum= np.zeros(readlen, np.int)
    nreads= np.zeros(readlen, np.int)

    for read in fq:
        if len(read)< len(qualsum):
            qualsum[:len(read)] += read.qual
            nreads[:len(read)] += 1
        else:
            qualsum += read.qual # it is a range of sum quality in each position [331 335 329 366 370 370 366 368 389 388 386 ...]
            nreads += 1
    print("number of reads :", nreads)
    print("quality sum  :", qualsum)

    y=qualsum/np.float_(nreads) # here we take the quality average [34 34 34 37 37 37 37 37 38 39 39 39 39 41 41 ...]
    print "average of quality in one file  :", y
    x = range(1,len(qualsum)+1,1)
    return(x,y)

def plotfigure( ftitle):
    plt.axis([0.0, 101, 20, 44])
    #ax = plt.gca()
    #ax.set_autoscale_on(False)
    plt.rc("font", size=16 , family ='serif')
    plt.title(ftitle, fontsize=26 , family ='serif')
    plt.ylabel('Phred score' ,fontsize=20 , family ='serif')
    plt.xlabel('Position in read' ,fontsize=20 , family ='serif')



samples = ['66604' , '66605' ,'66606' , '66607' ,'66608' , '66609' , '66610' , '66611' , '66613' ,'66614' , '66615' , '66616']
lebels = ['Day0_1' , 'Day0_2' , 'Day0_3' ,'Day3_1' ,'Day3_2' , 'Day3_3' ,'Day6_1' ,'Day6_2','Day6_3' ,'Day12_1' , 'Day12_2' , 'Day12_3']
colours = [ (0,0,0) , (0,1,0), (0 ,0, 1), (0,1,1 ) , (1,0,0), (1,0,1), (0.5,0,0.5), (0.5,0.5,0), (0.7,0,0.7), (0,0.4,0), (0.7,0.3,0), (0.9,0.6,0)]

# ------------------------------------------ After TRIMMING RNASeq ---------------------------------------------------------------------------------------------------------------------------------

fq_file_out = {}
sample_out = {}
fig2 = plt.figure(2, figsize=(8,8), dpi=120)
plotfigure('Mean Read Quality of QC_processed Data')

for i in range(len(samples)):
    print ("Sample %d: %s" % (i, samples[i]) )
    # dictionary of fastq sample
    sample_out["Forward" + str(i)]= pathout + samples[i] + '/' + samples[i] +'_1.fq.gz'
    sample_out["Reverse" + str(i)]= pathout + samples[i] + '/' + samples[i] +'_2.fq.gz'
    print(sample_out["Forward" + str(i)])
    print(sample_out["Reverse" + str(i)])
    #dictionary of FastqReader objects
    fq_file_out["Fq_F" + str(i)] = HTSeq.FastqReader ( sample_out["Forward" + str(i)])
    fq_file_out["Fq_R" + str(i)] = HTSeq.FastqReader ( sample_out["Reverse" + str(i)])
    # we passed the reverse and forward files to read_qual function to get the mean quality
    (Fxo,Fyo) = reads_qual(fq_file_out ["Fq_F" + str(i)])
    (Rxo,Ryo) = reads_qual(fq_file_out ["Fq_R" + str(i)])
    print "forward pos :" , Fxo
    print "Rev pos :" , Rxo
    # the average of forward mean quality and reverse mean quality
    xo = Fxo
    yo= (Fyo +Ryo)/float(2)
    print "average MQ :", yo
    print "length avrage quality", len(yo)
    plt.figure(2)
    plt.plot(xo ,yo ,color = colours[i],linestyle='solid',label= lebels[i], linewidth=2.0) # plot the mean quality of each sample
    plt.legend(loc='upper center', bbox_to_anchor=(0.56, 1.01),ncol=4, fancybox=True, shadow=True,fontsize=14 )

fig2.savefig( 'afterQC_quality' +'.png')
fig2.savefig( 'afterQC_quality'+'.pdf')