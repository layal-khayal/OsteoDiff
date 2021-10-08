__author__ = 'layal'

#this script calcules the coverage for the 12 samples by using the function of calculatecoverge

import HTSeq
import os
from calculateCoverage import calculateCoverageOfOneBamFile


# loop over all the BAM files all calculateCoverage
# change the path where is your data live
root= '/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO/Align_Tophat_GRCm38/'
gtf_files = '/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO/coverage/'

outCov_path=  '/media/layal/ac6ca8e5-d199-4d66-89a8-d1180b73dd9c/layal/OSTEO/coverage/coverage_output/'

if not os.path.exists(outCov_path):
    os.makedirs(outCov_path)


BAMname =[66604 , 66605 , 66606, 66607, 66608, 66609, 66610, 66611, 66613, 66614, 66615, 66616]
filesNameBase =['first_exon_cov_' , 'middle_exon_cov_', 'last_exon_cov_']

print "Read the gtf of the first exon"
gtf_file_firstExons = HTSeq.GFF_Reader(gtf_files+'first_exons_GRCm38.gtf')

print "Read the gtf of the last exon"
gtf_file_lastExons = HTSeq.GFF_Reader(gtf_files+'last_exons_GRCm38.gtf')

print "Read the gtf of the middle exon"
gtf_file_middleExons = HTSeq.GFF_Reader(gtf_files+'middle_exons_GRCm38.gtf')

#vector of all the gtf object
gtfFiles= [gtf_file_firstExons, gtf_file_middleExons, gtf_file_lastExons]
# for the first exons




for k, j in zip(gtfFiles , filesNameBase) :
    for i in BAMname:
        filename_cv =outCov_path + j +'%d.txt' %(i)
        print "creat alignment  object : "+'%d'%(i)
        Alignment = HTSeq.BAM_Reader(root+ str(i) + '/accepted_hits.bam')
        print "calculate the coverage from the annotation of " + j[:10] + " of the aligned sample :" +'%d'%(i)
        meanCov = calculateCoverageOfOneBamFile(k ,Alignment ,filename_cv)
        cov = open(outCov_path+'meanCov_'+j[:5]+str(i), 'w')
        cov.write(str(meanCov)+'\n')








