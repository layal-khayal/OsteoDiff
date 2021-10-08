# Author Layal
# version 2.0.0 last modified oct  2017
# Differential expression between RNA Seq Osteoblast time points
# Input data: BAM files mapped with tophat to the reference GRCm38 release 86
# annotation file is from ensemble GRCm38 release 86

#====================================================================
# INSTALL
# sudo apt-get install libcurl4-openssl-dev libxml2-dev
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# biocLite("IHW")
# install.packages("hash")
# biocLite("Rsubread")

library(Rsubread)
library(ggplot2)
library(DESeq2)
library(hash)

# to use parallel cores for faster calculation
library("BiocParallel")
register(MulticoreParam(8))

#the path to RNAseq data and tophat in SOLEXA
pathOsteo = "/home/layal/Solexa/layal/OSTEO"

#the output results path
pathOut = "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq"

#load(paste(path, "DESeq.RData", sep = "/"))

FC_countreads <- featureCounts(files=c(paste(pathOsteo, "Align_Tophat_GRCm38/66604/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66605/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66606/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66607/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66608/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66609/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66610/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66611/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66613/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66614/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66615/accepted_hits.bam", sep="/"),
                                       paste(pathOsteo, "Align_Tophat_GRCm38/66616/accepted_hits.bam", sep="/")),
                               isPairedEnd=TRUE, annot.ext= paste(pathOsteo,"GeneAnnotation/Mus_musculus.GRCm38.86.gtf",sep="/"),
                               isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_name")

saveRDS(FC_countreads , file = paste(pathOut,"FC_counts.rds",sep = "/"))
FC_countreads <- readRDS(paste(pathOut,"FC_counts.rds",sep = "/"))

# biomart table is downloaded from ensemble on 18
BioType <- read.table(paste(pathOut,"BioMart_GRCm38R86.txt", sep = "/") ,sep = "\t", header = T)

#FC_countreads_Ens <- readRDS(paste(path,"FC_counts_EnsemID.rds",sep = "/"))

countrawdata <- FC_countreads$counts
countdata <- countrawdata

colnames(countdata) <- gsub(pattern = "X.home.layal.Solexa.layal.OSTEO.Align_Tophat_GRCm38.|.accepted_hits.bam",replacement = "",x = colnames(countdata))

head(countdata)

New_names <- c("DayZero_R1","DayZero_R2","DayZero_R3" ,"DayThree_R1","DayThree_R2","DayThree_R3","DaySix_R1","DaySix_R2","DaySix_R3", "DayTwelve_R1", "DayTwelve_R2", "DayTwelve_R3")
colnames(countdata) <- New_names

#=================DEseq a nalysis : create DESeqDataSet ======================================================================================================================

condition <- factor(c(rep("DayZero", 3),rep("DayThree",3),rep("DaySix",3),rep("DayTwelve",3)), levels = c("DayZero", "DayThree","DaySix", "DayTwelve"))
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)

# **************=============================  DESeq function default ====================****************************
#removing rows in which there are no reads
ddsf <- dds[ rowSums(counts(dds)) > 1, ]
ddseq <- DESeq(ddsf)
assayNames(ddseq)

allgenes <- rownames(assay(ddseq))
write(allgenes, paste(pathOut,"allgenes.txt",sep = "/"))

#======================================== ** Normalization by Gene-specific normalization ** ==============================
#the gene-specific factors normalization, just genes specific
ddss<-ddsf
normFactors <- matrix(runif(nrow(ddss)*ncol(ddss),0.5,1.5),
                      ncol=ncol(ddss),nrow=nrow(ddss),
                      dimnames=list(1:nrow(ddss),1:ncol(ddss)))

normFactors <- normFactors / exp(rowMeans(log(normFactors)))

normalizationFactors(ddss) <- normFactors

# the subfunction of deseq() ::::::::
dSF <- estimateSizeFactors(ddss)
DisNF <- estimateDispersions(dSF)
desqNF <- nbinomWaldTest(DisNF)

assayNames(desqNF)
sizeFactors(desqNF)
normalizationFactors(desqNF)[1:3,]

#============================================================================== counts ===============================================================================

RawCounts <- counts(ddseq, normalized=FALSE)
write.csv(RawCounts, paste(pathOut,"rawCounts.csv",sep = "/"))

NormSF_Counts <- counts(ddseq, normalized=TRUE)
write.csv(NormSF_Counts, paste(pathOut,"NormSF_Counts.csv",sep = "/"))

NormSNF_Counts <- counts(desqNF, normalized=TRUE)
write.csv(NormSF_Counts, paste(pathOut,"NormSNF_Counts.csv",sep = "/"))

#=========================================================== Results ======================================================

# Day 12 to Day 0 comparison :

resTwelve_Zero <- results(ddseq, contrast=c("condition", "DayTwelve", "DayZero"),pAdjustMethod = "bonferroni", alpha = 0.01,
                          cooksCutoff = TRUE, independentFiltering = TRUE)
summary(resTwelve_Zero)
res12_0 <- as.data.frame(resTwelve_Zero)
for (i in 1:length(rownames(res12_0 ))) {
  res12_0$Gene.type[i] <- as.character(hash_biotype[[ rownames(res12_0)[i] ]])
}

resSix_Zero <- results(ddseq, contrast=c("condition", "DaySix" ,"DayZero") , pAdjustMethod = "bonferroni", alpha = 0.01,
                       cooksCutoff = TRUE, independentFiltering = TRUE)
res6_0 <- as.data.frame(resSix_Zero)
for (i in 1:length(rownames(res6_0 ))) {
  res6_0$Gene.type[i] <- as.character(hash_biotype[[ rownames(res6_0)[i] ]])
}

resThree_Zero <- results(ddseq, contrast=c("condition", "DayThree", "DayZero"), pAdjustMethod = "bonferroni", alpha = 0.01,
                         cooksCutoff = TRUE, independentFiltering = TRUE)
res3_0 <- as.data.frame(resThree_Zero)
for (i in 1:length(rownames(res3_0 ))) {
  res3_0$Gene.type[i] <- as.character(hash_biotype[[ rownames(res3_0)[i] ]])
}

resTwelve_Three <- results(ddseq, contrast=c("condition", "DayTwelve", "DayThree"),pAdjustMethod = "bonferroni", alpha = 0.01,
                           cooksCutoff = TRUE, independentFiltering = TRUE)
res12_3 <- as.data.frame(resTwelve_Three)
for (i in 1:length(rownames(res12_3 ))) {
  res12_3$Gene.type[i] <- as.character(hash_biotype[[ rownames(res12_3)[i] ]])
}

resTwelve_Six <- results(ddseq, contrast=c("condition", "DayTwelve", "DaySix"),pAdjustMethod = "bonferroni", alpha = 0.01,
                         cooksCutoff = TRUE, independentFiltering = TRUE)
res12_6 <- as.data.frame(resTwelve_Six)
for (i in 1:length(rownames(res12_6 ))) {
  res12_6$Gene.type[i] <- as.character(hash_biotype[[ rownames(res12_6)[i] ]])
}

resSix_Three <- results(ddseq, contrast=c("condition", "DaySix", "DayThree"),pAdjustMethod = "bonferroni", alpha = 0.01,
                        cooksCutoff = TRUE, independentFiltering = TRUE)
res6_3 <- as.data.frame(resSix_Three)
for (i in 1:length(rownames(res6_3 ))) {
  res6_3$Gene.type[i] <- as.character(hash_biotype[[ rownames(res6_3)[i] ]])
}
#=============================================================DEFAULT GENE EXPRESION THRESHOLD ======================
pd = 0.01
lfc = log2(4)
#-------------------------12_0-----------------------------------
#10058
sig12_0 <- res12_0 [ which(res12_0$padj< pd),]
#2299
Diff12_0 <- res12_0 [ which(res12_0$padj < pd & abs(res12_0$log2FoldChange) > lfc),]
#1231
Diff12_0_up <- res12_0 [ which(res12_0$padj < pd & res12_0$log2FoldChange > lfc),]
#1068
Diff12_0_dn <- res12_0 [ which(res12_0$padj < pd & res12_0$log2FoldChange < -lfc),]

#-------------------------6_0-----------------------------------
#10026
sig6_0 <- res6_0 [ which(res6_0$padj< pd),]
#2834
Diff6_0 <- res6_0 [ which(res6_0$padj < pd & abs(res6_0$log2FoldChange) > lfc),]
#1702
Diff6_0_up <- res6_0 [ which(res6_0$padj < pd & res6_0$log2FoldChange > lfc),]
#1132
Diff6_0_dn <- res6_0 [ which(res6_0$padj < pd & res6_0$log2FoldChange < -lfc),]

#-------------------------3_0-----------------------------------
#9695
sig3_0 <- res3_0 [ which(res3_0$padj< pd),]
#2378
Diff3_0 <- res3_0 [ which(res3_0$padj < pd & abs(res3_0$log2FoldChange) > lfc),]
#1466
Diff3_0_up <- res3_0 [ which(res3_0$padj < pd & res3_0$log2FoldChange > lfc),]
#912
Diff3_0_dn <- res3_0 [ which(res3_0$padj < pd & res3_0$log2FoldChange < -lfc),]

#-------------------------12_3-----------------------------------
#5551
sig12_3 <- res12_3 [ which(res12_3$padj< pd),]
#655
Diff12_3 <- res12_3 [ which(res12_3$padj < pd & abs(res12_3$log2FoldChange) > lfc),]
#233
Diff12_3_up <- res12_3 [ which(res12_3$padj < pd & res12_3$log2FoldChange > lfc),]
#422
Diff12_3_dn <- res12_3 [ which(res12_3$padj < pd & res12_3$log2FoldChange < -lfc),]

#-------------------------12_6-----------------------------------
#5427
sig12_6 <- res12_6 [ which(res12_6$padj< pd),]
#626
Diff12_6 <- res12_6 [ which(res12_6$padj < pd & abs(res12_6$log2FoldChange) > lfc),]
#214
Diff12_6_up <- res12_6 [ which(res12_6$padj < pd & res12_6$log2FoldChange > lfc),]
#412
Diff12_6_dn <- res12_6 [ which(res12_6$padj < pd & res12_6$log2FoldChange < -lfc),]

#-------------------------6_3-----------------------------------
#1788
sig6_3 <-  res6_3 [ which(res6_3$padj < pd),]
#155
Diff6_3 <- res6_3 [ which(res6_3$padj < pd & abs(res6_3$log2FoldChange) > lfc),]
#92
Diff6_3up <- res6_3 [ which(res6_3$padj < pd & res6_3$log2FoldChange > lfc),]
#63
Diff6_3_dn <- res6_3 [ which(res6_3$padj < pd & res6_3$log2FoldChange < -lfc),]

#==========================================================================biotype====================================================================================================
myBioType <- subset(BioType, BioType$Associated.Gene.Name %in% rownames(ddseq))
#without the ensemble ID
pureBioType <- myBioType[,-1]
#detete the dublicated genes
pureBioType1 <- pureBioType[-which(duplicated(pureBioType$Associated.Gene.Name) == TRUE),]

#----------preparing the biotype ---------

## detete the dublicated genes, but some of them are protein_coding but we delete them and took the second biotype so we need tocheck this
# first find what are the gene which are duplicated, they are 158 genes
test1 <- pureBioType[which(duplicated(myBioType$Associated.Gene.Name) == TRUE),]
#second take the subset of all of them
test2 <- subset(pureBioType, pureBioType$Associated.Gene.Name %in% test1$Associated.Gene.Name)
#third the duplicated gene with protein_coding
tt<-test2[which(test2$Gene.type == "protein_coding"),]
td <-tt[which(duplicated(tt$Associated.Gene.Name) == TRUE),]
ts <- subset(tt, tt$Associated.Gene.Name %in% td$Associated.Gene.Name)
# substract the duplicated gene name with protein_ coding biotype from tt , we will find the genes their name are duplicated and have protein coding biotype and another one
tm <- tt[! tt$Associated.Gene.Name %in% ts$Associated.Gene.Name , ]
#the gene we must check that there biotype is protein coding
genestomod <- as.character(tm$Associated.Gene.Name)
pureBioType2 <- pureBioType1
for (i in genestomod) {
  pureBioType2[which(pureBioType2$Associated.Gene.Name == i ),2]<-"protein_coding"
}

#------check the number of protien coding -------------
#17945
p1 <- pureBioType1[which(pureBioType1$Gene.type == "protein_coding"),]
#17948
p2 <- pureBioType2[which(pureBioType2$Gene.type == "protein_coding"),]
#11200
np <- pureBioType2[which(pureBioType2$Gene.type != "protein_coding"),]
#2008
lncR <- pureBioType2[which(pureBioType2$Gene.type == "lincRNA"),]

hash_biotype <- hash(keys= pureBioType2$Associated.Gene.Name , values= pureBioType2$Gene.type) 

#======================================= THE THRESHOLD FOR THE HEATMAP==============================================
hash_Pvalue3_0 <- hash(keys= rownames(res3_0) , values= res3_0$padj)
hash_Pvalue6_0 <- hash(keys= rownames(res6_0) , values= res6_0$padj)
hash_Pvalue12_0 <- hash(keys= rownames(res12_0) , values= res12_0$padj)
hash_Pvalue6_3 <- hash(keys= rownames(res6_3) , values= res6_3$padj)
hash_Pvalue12_3 <- hash(keys= rownames(res12_3) , values= res12_3$padj)
hash_Pvalue12_6 <- hash(keys= rownames(res12_6) , values= res12_6$padj)

hash_LFC3_0 <- hash(keys= rownames(res3_0) , values= res3_0$log2FoldChange)
hash_LFC6_0 <- hash(keys= rownames(res6_0) , values= res6_0$log2FoldChange)
hash_LFC12_0 <- hash(keys= rownames(res12_0) , values= res12_0$log2FoldChange)
hash_LFC6_3 <- hash(keys= rownames(res6_3) , values= res6_3$log2FoldChange)
hash_LFC12_3 <- hash(keys= rownames(res12_3) , values= res12_3$log2FoldChange)
hash_LFC12_6 <- hash(keys= rownames(res12_6) , values= res12_6$log2FoldChange)
 
DifG <- data.frame(row.names = rownames(ddseq))
DifGA <- DifG
for (i in 1:length(rownames(DifGA))) {
  DifGA$LFC3_0[i] <- hash_LFC3_0[[ rownames(DifGA)[i] ]]
  DifGA$LFC6_0[i] <- hash_LFC6_0[[ rownames(DifGA)[i] ]]
  DifGA$LFC12_0[i] <- hash_LFC12_0[[ rownames(DifGA)[i] ]]
  DifGA$LFC12_3[i] <- hash_LFC12_3[[ rownames(DifGA)[i] ]]
  DifGA$LFC12_6[i] <- hash_LFC12_6[[ rownames(DifGA)[i] ]]
  DifGA$LFC6_3[i] <- hash_LFC6_3[[ rownames(DifGA)[i] ]]
  
  DifGA$Padj3_0[i] <- hash_Pvalue3_0[[ rownames(DifGA)[i] ]]
  DifGA$Padj6_0[i] <- hash_Pvalue6_0[[ rownames(DifGA)[i] ]]
  DifGA$Padj12_0[i] <- hash_Pvalue12_0[[ rownames(DifGA)[i] ]]
  DifGA$Padj12_3[i] <- hash_Pvalue12_3[[ rownames(DifGA)[i] ]]
  DifGA$Padj12_6[i] <- hash_Pvalue12_6[[ rownames(DifGA)[i] ]]
  DifGA$Padj6_3[i] <- hash_Pvalue6_3[[ rownames(DifGA)[i] ]]
  
  DifGA$Gene.type[i] <- as.character(hash_biotype[[ rownames(DifGA)[i] ]])
}


#==============================================================p-value < E-50============================================
my_FC <- log2(5)
my_padj <- 10^(-50)

#4537
SigE50 <- DifGA[ which( (DifGA$Padj3_0 < my_padj  ) | (DifGA$Padj6_0 < my_padj) | (DifGA$Padj12_0 < my_padj ) | (DifGA$Padj12_3 < my_padj )
                        | (DifGA$Padj12_6 < my_padj ) | (DifGA$Padj6_3 < my_padj)), ]
#4449
sigPro_E50 <- SigE50[which(SigE50$Gene.type == "protein_coding"),]
#88
sigNonCod_E50 <- SigE50[which(SigE50$Gene.type != "protein_coding"),]
#33
siglncRNA_E50 <- SigE50[which(SigE50$Gene.type == "lincRNA"),]


#1441
DiffGenesAll_E_50 <- DifGA[ which( (DifGA$Padj3_0 < my_padj & abs(DifGA$LFC3_0) > my_FC ) | (DifGA$Padj6_0 < my_padj & abs(DifGA$LFC6_0) > my_FC )
                                   | (DifGA$Padj12_0 < my_padj & abs(DifGA$LFC12_0) > my_FC ) | (DifGA$Padj12_3 < my_padj & abs(DifGA$LFC12_3) > my_FC )
                                   | (DifGA$Padj12_6 < my_padj & abs(DifGA$LFC12_6) > my_FC ) | (DifGA$Padj6_3 < my_padj & abs(DifGA$LFC6_3) > my_FC )), ]
#826
DiffGenesAll_E_50_up <- DifGA[ which( (DifGA$Padj3_0 < my_padj & DifGA$LFC3_0 > my_FC ) | (DifGA$Padj6_0 < my_padj & DifGA$LFC6_0 > my_FC )
                                   | (DifGA$Padj12_0 < my_padj & DifGA$LFC12_0 > my_FC ) | (DifGA$Padj12_3 < my_padj & DifGA$LFC12_3 > my_FC )
                                   | (DifGA$Padj12_6 < my_padj & DifGA$LFC12_6 > my_FC ) | (DifGA$Padj6_3 < my_padj & DifGA$LFC6_3 > my_FC )), ]
#781
DiffGenesAll_E_50_dn <- DifGA[ which( (DifGA$Padj3_0 < my_padj & DifGA$LFC3_0 < -my_FC ) | (DifGA$Padj6_0 < my_padj & DifGA$LFC6_0 < -my_FC )
                                      | (DifGA$Padj12_0 < my_padj & DifGA$LFC12_0 < -my_FC ) | (DifGA$Padj12_3 < my_padj & DifGA$LFC12_3 < -my_FC )
                                      | (DifGA$Padj12_6 < my_padj & DifGA$LFC12_6 < -my_FC ) | (DifGA$Padj6_3 < my_padj & DifGA$LFC6_3 < -my_FC )), ]

#1386 (the genes for the heatmap)
DiffGenesPro_E50 <- DiffGenesAll_E_50[which(DiffGenesAll_E_50$Gene.type == "protein_coding"),]
#55
DiffGenesNonCod_E50 <- DiffGenesAll_E_50[which(DiffGenesAll_E_50$Gene.type != "protein_coding"),]
#19
DiffGenes_lncRNA_E50 <- DiffGenesAll_E_50[which(DiffGenesAll_E_50$Gene.type == "lincRNA"),]

#===============================================================================================
pv <- 0.01
fc <- log2(4)

sig0.01 <- DifGA[ which( (DifGA$Padj3_0 < pv) | (DifGA$Padj6_0 < pv ) | (DifGA$Padj12_0 < pv ) | (DifGA$Padj12_3 < pv )
                      | (DifGA$Padj12_6 < pv ) | (DifGA$Padj6_3 < pv )), ]
#11773
sigPro0.01<- sig0.01[which(sig0.01$Gene.type == "protein_coding"),]
#1159
sigNonCod_0.01 <- sig0.01[which(sig0.01$Gene.type != "protein_coding"),]
#287
siglncRNA_0.01 <- sig0.01[which(sig0.01$Gene.type == "lincRNA"),]


#4012
DifGE <- DifGA[ which( (DifGA$Padj3_0 < pv & abs(DifGA$LFC3_0) > fc ) | (DifGA$Padj6_0 < pv & abs(DifGA$LFC6_0) > fc )
                    | (DifGA$Padj12_0 < pv & abs(DifGA$LFC12_0) > fc ) | (DifGA$Padj12_3 < pv & abs(DifGA$LFC12_3) > fc )
                    | (DifGA$Padj12_6 < pv & abs(DifGA$LFC12_6) > fc ) | (DifGA$Padj6_3 < pv & abs(DifGA$LFC6_3) > fc )), ]

#2414
DifGE_up <- DifGA[ which( (DifGA$Padj3_0 < pv & DifGA$LFC3_0 > fc ) | (DifGA$Padj6_0 < pv & DifGA$LFC6_0 > fc )
                       | (DifGA$Padj12_0 < pv & DifGA$LFC12_0 > fc ) | (DifGA$Padj12_3 < pv & DifGA$LFC12_3 > fc )
                       | (DifGA$Padj12_6 < pv & DifGA$LFC12_6 > fc ) | (DifGA$Padj6_3 < pv & DifGA$LFC6_3 > fc )), ]

#2124
DifGE_dn <- DifGA[ which( (DifGA$Padj3_0 < pv & DifGA$LFC3_0 < -fc ) | (DifGA$Padj6_0 < pv & DifGA$LFC6_0 < -fc )
                          | (DifGA$Padj12_0 < pv & DifGA$LFC12_0 < -fc ) | (DifGA$Padj12_3 < pv & DifGA$LFC12_3 < -fc )
                          | (DifGA$Padj12_6 < pv & DifGA$LFC12_6 < -fc ) | (DifGA$Padj6_3 < pv & DifGA$LFC6_3 < -fc )), ]


#3308
DifGE_pro <- DifGE[which(DifGE$Gene.type == "protein_coding"),]
#704
DifGE_nonCod <- DifGE[which(DifGE$Gene.type != "protein_coding"),]
#184
DifGE_linc <- DifGE[which(DifGE$Gene.type == "lincRNA"),]

#================================================ linc===============================================================
pvlnc <- 0.01
fclnc <- log2(1.5)

#11928
Diflnc <- DifGA[ which( (DifGA$Padj3_0 < pvlnc & abs(DifGA$LFC3_0) > fclnc ) | (DifGA$Padj6_0 < pvlnc & abs(DifGA$LFC6_0) > fclnc )
                       | (DifGA$Padj12_0 < pvlnc & abs(DifGA$LFC12_0) > fclnc ) | (DifGA$Padj12_3 < pvlnc & abs(DifGA$LFC12_3) > fclnc )
                       | (DifGA$Padj12_6 < pvlnc & abs(DifGA$LFC12_6) > fclnc ) | (DifGA$Padj6_3 < pvlnc & abs(DifGA$LFC6_3) > fclnc )), ]

#10778
Diflnc_pro <- Diflnc[which(Diflnc$Gene.type == "protein_coding"),]
#1150
Diflnc_nonCod <- Diflnc[which(Diflnc$Gene.type != "protein_coding"),]
#285
Diflnc_linc <- Diflnc[which(Diflnc$Gene.type == "lincRNA"),]

pathTables = "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq/DESeq_table"
# exporting table
write.csv(DifGE,paste(pathTables,"DGEFC4P0.01.csv",sep = "/"))
write.csv(DifGE_pro,paste(pathTables,"DGE_pro.csv",sep = "/"))
write.csv(DifGE_nonCod,paste(pathTables,"DGE_nonCod.csv",sep = "/"))

# tables of comparisons DGE up and down
#------------------12_0-------------------------------------

write.csv(Diff12_0_up,paste(pathTables,"Diff12_0Up.csv",sep = "/"))
write.csv(Diff12_0_dn,paste(pathTables,"Diff12_0Dn.csv",sep = "/"))

#----------------------6_0-------------------------------

write.csv(Diff6_0_up,paste(pathTables,"Diff6_0Up.csv",sep = "/"))
write.csv(Diff6_0_dn,paste(pathTables,"Diff6_0Dn.csv",sep = "/"))

#----------------------3_0-------------------------------

write.csv(Diff3_0_up,paste(pathTables,"Diff3_0Up.csv",sep = "/"))
write.csv(Diff3_0_dn,paste(pathTables,"Diff3_0Dn.csv",sep = "/"))
#-----------------------12_3-----------------

write.csv(Diff12_3_up,paste(pathTables,"Diff12_3Up.csv",sep = "/"))
write.csv(Diff12_3_dn,paste(pathTables,"Diff12_3Dn.csv",sep = "/"))

#----------------------12_6-------------------------------
write.csv(Diff12_6_up,paste(pathTables,"Diff12_6Up.csv",sep = "/"))
write.csv(Diff12_6_dn,paste(pathTables,"Diff12_6Dn.csv",sep = "/"))
#--------------------6_3----------------------------------------------
write.csv(Diff6_3up,paste(pathTables,"Diff6_3Up.csv",sep = "/"))
write.csv(Diff6_3_dn,paste(pathTables,"Diff6_3Dn.csv",sep = "/"))

# in this format R read the file correctly as matrix by read.table
write.table(DiffGenesPro_E50 ,paste(pathTables,"DiffGenesPro_E50.txt",sep = "/"), sep = "\t")
write.csv(DiffGenesPro_E50 ,paste(pathTables,"DiffGenesPro_E50.csv",sep = "/"))
write.table(DiffGenesNonCod_E50 ,paste(pathTables,"DiffGenesNonCod_E50.txt",sep = "/"), sep = "\t")
write.csv(DiffGenesNonCod_E50 ,paste(pathTables,"DiffGenesNoncod_E50.csv",sep = "/"))

genesProt <- rownames(DiffGenesPro_E50)
write(genesProt, paste(pathTables,"genesProt.txt",sep = "/"))

genesNonCod <- rownames(DiffGenesNonCod_E50)
write(genesNonCod, paste(pathTables,"genesNonCod.txt",sep = "/"))

#================================ check the on off groups between samples ============================
ON_Off <- data.frame()
#---------------------------------------- Day 3_0-------------------------------------------------------------------

NormCoun3_0 <- as.data.frame( NormSF_Counts[which(rownames(NormSF_Counts) %in% rownames(Diff3_0)), 1:6])

NormCoun3_0$Mean0 <- rowMeans(cbind(NormCoun3_0$DayZero_R1,NormCoun3_0$DayZero_R2,NormCoun3_0$DayZero_R3))
NormCoun3_0$Mean3 <- rowMeans(cbind(NormCoun3_0$DayThree_R1,NormCoun3_0$DayThree_R2,NormCoun3_0$DayThree_R3))

for (i in 1:length(rownames(NormCoun3_0))) {
  NormCoun3_0$Gene.type[i] <- as.character(DifBiotyp_hash[[ rownames(NormCoun3_0)[i] ]])
}


On3off0 <-  NormCoun3_0[which((NormCoun3_0$DayZero_R1 ==0 &NormCoun3_0$DayZero_R2 == 0 & NormCoun3_0$Mean0 < 1.5) |
                                (NormCoun3_0$DayZero_R1 ==0 &NormCoun3_0$DayZero_R3 == 0 & NormCoun3_0$Mean0 < 1.5) |
                                (NormCoun3_0$DayZero_R2 ==0 &NormCoun3_0$DayZero_R3 == 0 & NormCoun3_0$Mean0 < 1.5) ), ]

On0off3 <-  NormCoun3_0[which((NormCoun3_0$DayThree_R1 ==0 &NormCoun3_0$DayThree_R2 == 0 & NormCoun3_0$Mean3 < 1.5) |
                                (NormCoun3_0$DayThree_R1 ==0 &NormCoun3_0$DayThree_R3 == 0 & NormCoun3_0$Mean3 < 1.5) |
                                (NormCoun3_0$DayThree_R2 ==0 &NormCoun3_0$DayThree_R3 == 0 & NormCoun3_0$Mean3 < 1.5) ), ]

write.csv(On3off0,paste(pathTables,"on3off0.csv",sep = "/"))
write.csv(On0off3,paste(pathTables,"on0off3.csv",sep = "/"))


#------------------------------------------- Day 6_0 --------------------------------------------------------------------------

NormCoun6_0 <- as.data.frame( NormSF_Counts[which(rownames(NormSF_Counts) %in% rownames(Diff6_0)), c(1:3,7:9)])

NormCoun6_0$Mean0 <- rowMeans(cbind(NormCoun6_0$DayZero_R1,NormCoun6_0$DayZero_R2,NormCoun6_0$DayZero_R3))
NormCoun6_0$Mean6 <- rowMeans(cbind(NormCoun6_0$DaySix_R1,NormCoun6_0$DaySix_R2,NormCoun6_0$DaySix_R3))

for (i in 1:length(rownames(NormCoun6_0))) {
  NormCoun6_0$Gene.type[i] <- as.character(DifBiotyp_hash[[ rownames(NormCoun6_0)[i] ]])
}

On6off0 <-  NormCoun6_0[which((NormCoun6_0$DayZero_R1 ==0 &NormCoun6_0$DayZero_R2 == 0 & NormCoun6_0$Mean0 < 1.5) |
                                (NormCoun6_0$DayZero_R1 ==0 &NormCoun6_0$DayZero_R3 == 0 & NormCoun6_0$Mean0 < 1.5) |
                                (NormCoun6_0$DayZero_R2 ==0 &NormCoun6_0$DayZero_R3 == 0 & NormCoun6_0$Mean0 < 1.5) ), ]

On0off6 <-  NormCoun6_0[which((NormCoun6_0$DaySix_R1 ==0 &NormCoun6_0$DaySix_R2 == 0 & NormCoun6_0$Mean6 < 1.5) |
                                (NormCoun6_0$DaySix_R1 ==0 &NormCoun6_0$DaySix_R3 == 0 & NormCoun6_0$Mean6 < 1.5) |
                                (NormCoun6_0$DaySix_R2 ==0 &NormCoun6_0$DaySix_R3 == 0 & NormCoun6_0$Mean6 < 1.5) ), ]

write.csv(On6off0,paste(pathTables,"on6off0.csv",sep = "/"))
write.csv(On0off6,paste(pathTables,"on0off6.csv",sep = "/"))

#--------------------------------------------------- Day 12_0 --------------------------------------------------

NormCoun12_0 <- as.data.frame( NormSF_Counts[which(rownames(NormSF_Counts) %in% rownames(Diff12_0)), c(1:3,10:12)])

NormCoun12_0$Mean0 <- rowMeans(cbind(NormCoun12_0$DayZero_R1,NormCoun12_0$DayZero_R2,NormCoun12_0$DayZero_R3))
NormCoun12_0$Mean12 <- rowMeans(cbind(NormCoun12_0$DayTwelve_R1,NormCoun12_0$DayTwelve_R2,NormCoun12_0$DayTwelve_R3))

for (i in 1:length(rownames(NormCoun12_0))) {
  NormCoun12_0$Gene.type[i] <- as.character(DifBiotyp_hash[[ rownames(NormCoun12_0)[i] ]])
}


On12off0 <-  NormCoun12_0[which((NormCoun12_0$DayZero_R1 ==0 &NormCoun12_0$DayZero_R2 == 0 & NormCoun12_0$Mean0 < 1.5) |
                                (NormCoun12_0$DayZero_R1 ==0 &NormCoun12_0$DayZero_R3 == 0 & NormCoun12_0$Mean0 < 1.5) |
                                (NormCoun12_0$DayZero_R2 ==0 &NormCoun12_0$DayZero_R3 == 0 & NormCoun12_0$Mean0 < 1.5) ), ]

On0off12 <-  NormCoun12_0[which((NormCoun12_0$DayTwelve_R1 ==0 &NormCoun12_0$DayTwelve_R2 == 0 & NormCoun12_0$Mean12 < 1.5) |
                                (NormCoun12_0$DayTwelve_R1 ==0 &NormCoun12_0$DayTwelve_R3 == 0 & NormCoun12_0$Mean12 < 1.5) |
                                (NormCoun12_0$DayTwelve_R2 ==0 &NormCoun12_0$DayTwelve_R3 == 0 & NormCoun12_0$Mean12 < 1.5) ), ]

write.csv(On12off0,paste(pathTables,"on12off0.csv",sep = "/"))
write.csv(On0off12,paste(pathTables,"on0off12.csv",sep = "/"))

#------------------------------------------------------ Day 12_3----------------------------------------------------

NormCoun12_3 <- as.data.frame( NormSF_Counts[which(rownames(NormSF_Counts) %in% rownames(Diff12_3)), c(4:6,10:12)])

NormCoun12_3$Mean3 <- rowMeans(cbind(NormCoun12_3$DayThree_R1,NormCoun12_3$DayThree_R2,NormCoun12_3$DayThree_R3))
NormCoun12_3$Mean12 <- rowMeans(cbind(NormCoun12_3$DayTwelve_R1,NormCoun12_3$DayTwelve_R2,NormCoun12_3$DayTwelve_R3))

for (i in 1:length(rownames(NormCoun12_3))) {
  NormCoun12_3$Gene.type[i] <- as.character(DifBiotyp_hash[[ rownames(NormCoun12_3)[i] ]])
}

On12off3 <-  NormCoun12_3[which((NormCoun12_3$DayThree_R1 ==0 &NormCoun12_3$DayThree_R2 == 0 & NormCoun12_3$Mean3 < 1.5) |
                                  (NormCoun12_3$DayThree_R1 ==0 &NormCoun12_3$DayThree_R3 == 0 & NormCoun12_3$Mean3 < 1.5) |
                                  (NormCoun12_3$DayThree_R2 ==0 &NormCoun12_3$DayThree_R3 == 0 & NormCoun12_3$Mean3 < 1.5) ), ]


On3off12 <-  NormCoun12_3[which((NormCoun12_3$DayTwelve_R1 ==0 &NormCoun12_3$DayTwelve_R2 == 0 & NormCoun12_3$Mean12 < 1.5) |
                                  (NormCoun12_3$DayTwelve_R1 ==0 &NormCoun12_3$DayTwelve_R3 == 0 & NormCoun12_3$Mean12 < 1.5) |
                                  (NormCoun12_3$DayTwelve_R2 ==0 &NormCoun12_3$DayTwelve_R3 == 0 & NormCoun12_3$Mean12 < 1.5) ), ]

write.csv(On12off3,paste(pathTables,"on12off3.csv",sep = "/"))
write.csv(On3off12,paste(pathTables,"on3off12.csv",sep = "/"))

#--------------------------------------------------- Day 12_6 ------------------------------------------------

NormCoun12_6 <- as.data.frame( NormSF_Counts[which(rownames(NormSF_Counts) %in% rownames(Diff12_6)), 7:12])

NormCoun12_6$Mean6 <- rowMeans(cbind(NormCoun12_6$DaySix_R1,NormCoun12_6$DaySix_R2,NormCoun12_6$DaySix_R3))
NormCoun12_6$Mean12 <- rowMeans(cbind(NormCoun12_6$DayTwelve_R1,NormCoun12_6$DayTwelve_R2,NormCoun12_6$DayTwelve_R3))

for (i in 1:length(rownames(NormCoun12_6))) {
  NormCoun12_6$Gene.type[i] <- as.character(DifBiotyp_hash[[ rownames(NormCoun12_6)[i] ]])
}


On12off6 <-  NormCoun12_6[which((NormCoun12_6$DaySix_R1 ==0 &NormCoun12_6$DaySix_R2 == 0 & NormCoun12_6$Mean6 < 1.5) |
                                  (NormCoun12_6$DaySix_R1 ==0 &NormCoun12_6$DaySix_R3 == 0 & NormCoun12_6$Mean6 < 1.5) |
                                  (NormCoun12_6$DaySix_R2 ==0 &NormCoun12_6$DaySix_R3 == 0 & NormCoun12_6$Mean6 < 1.5) ), ]


On6off12 <-  NormCoun12_6[which((NormCoun12_6$DayTwelve_R1 ==0 &NormCoun12_6$DayTwelve_R2 == 0 & NormCoun12_6$Mean12 < 1.5) |
                                  (NormCoun12_6$DayTwelve_R1 ==0 &NormCoun12_6$DayTwelve_R3 == 0 & NormCoun12_6$Mean12 < 1.5) |
                                  (NormCoun12_6$DayTwelve_R2 ==0 &NormCoun12_6$DayTwelve_R3 == 0 & NormCoun12_6$Mean12 < 1.5) ), ]

write.csv(On12off6,paste(pathTables,"on12off6.csv",sep = "/"))
write.csv(On6off12,paste(pathTables,"on6off12.csv",sep = "/"))

#------------------------------------------------------ Day 6_3----------------------------------------------------

NormCoun6_3 <- as.data.frame( NormSF_Counts[which(rownames(NormSF_Counts) %in% rownames(Diff6_3)), 4:9])

NormCoun6_3$Mean3 <- rowMeans(cbind(NormCoun6_3$DayThree_R1,NormCoun6_3$DayThree_R2,NormCoun6_3$DayThree_R3))
NormCoun6_3$Mean6 <- rowMeans(cbind(NormCoun6_3$DaySix_R1,NormCoun6_3$DaySix_R2,NormCoun6_3$DaySix_R3))

for (i in 1:length(rownames(NormCoun6_3))) {
  NormCoun6_3$Gene.type[i] <- as.character(DifBiotyp_hash[[ rownames(NormCoun6_3)[i] ]])
}


On6off3 <-  NormCoun6_3[which((NormCoun6_3$DayThree_R1 ==0 &NormCoun6_3$DayThree_R2 == 0 & NormCoun6_3$Mean3 < 1.5) |
                                  (NormCoun6_3$DayThree_R1 ==0 &NormCoun6_3$DayThree_R3 == 0 & NormCoun6_3$Mean3 < 1.5) |
                                  (NormCoun6_3$DayThree_R2 ==0 &NormCoun6_3$DayThree_R3 == 0 & NormCoun6_3$Mean3 < 1.5) ), ]

On3off6 <-  NormCoun6_3[which((NormCoun6_3$DaySix_R1 ==0 &NormCoun6_3$DaySix_R2 == 0 & NormCoun6_3$Mean6 < 1.5) |
                                  (NormCoun6_3$DaySix_R1 ==0 &NormCoun6_3$DaySix_R3 == 0 & NormCoun6_3$Mean6 < 1.5) |
                                  (NormCoun6_3$DaySix_R2 ==0 &NormCoun6_3$DaySix_R3 == 0 & NormCoun6_3$Mean6 < 1.5) ), ]

write.csv(On6off3,paste(pathTables,"on6off3.csv",sep = "/"))
write.csv(On3off6,paste(pathTables,"on3off6.csv",sep = "/"))

# ----------------------- table on off union-----------------------------------------------------
OnOffGenes <- union(rownames(On0off3), union(rownames(On0off6),union( rownames(On0off12),
                     union(rownames(On3off0),union(rownames(On3off6),union(rownames(On3off12),
                    union( rownames(On6off12), union(rownames(On6off3), union(rownames(On6off0),
                    union(rownames(On12off6), union(rownames(On12off0), rownames(On12off3))))))))))))

annot <- read.table("/home/layal/OSTEO_2017/DESeq/biomart86TopFun.txt" ,sep = "\t", header = T)
ant <- annot[which(annot$Associated.Gene.Name %in% OnOffGenes), ]
hash_biotype <- hash(keys= ant$Associated.Gene.Name , values= ant$Gene.type)
hash_G.start <- hash(keys= ant$Associated.Gene.Name , values= ant$Gene.Start..bp.)
hash_G.end   <- hash(keys= ant$Associated.Gene.Name , values= ant$Gene.End..bp.)
hash_chr <- hash(keys= ant$Associated.Gene.Name, values=ant$Chromosome.Name)


NCOF <- as.data.frame(NormSF_Counts[which(rownames(NormSF_Counts) %in% OnOffGenes), ]) 

NCOF$meanZ <- rowMeans(cbind(NCOF$DayZero_R1, NCOF$DayZero_R2 , NCOF$DayZero_R3), na.rm=TRUE)
NCOF$ meanTh <- rowMeans(cbind(NCOF$DayThree_R1 , NCOF$DayThree_R2 , NCOF$DayThree_R3) , na.rm = TRUE)
NCOF$meanS <- rowMeans(cbind(NCOF$DaySix_R1 , NCOF$DaySix_R2, NCOF$DaySix_R3) , na.rm = TRUE)
NCOF$meanTw <- rowMeans(cbind(NCOF$DayTwelve_R1, NCOF$DayTwelve_R2 , NCOF$DayTwelve_R3) , na.rm = TRUE)

COFG <- NCOF[,13:16]

for(i in 1:length(rownames(COFG))){
  COFG$Biotype[i] <- as.character(hash_biotype[[rownames(COFG)[i] ]])
  COFG$chr[i] <- as.character(hash_chr[[rownames(COFG)[i] ]])
  COFG$start[i] <- hash_G.start[[rownames(COFG)[i] ]]
  COFG$end[i] <- hash_G.end[[ rownames(COFG)[i] ]]
}

n <- length(rownames(COFG))
COFG$Mean_0 <- rep(" ", n)
COFG$Mean_3 <- rep(" ", n)
COFG$Mean_6 <- rep(" ", n)
COFG$Mean_12 <- rep(" ", n)

COFG <- COFG[,c(5:8,1:4,9:12)] 

for(i in 1:length(rownames(COFG))){
  L <- COFG[i,5:8]
  pk <- max(L)
  low <- which(L > 1.5 & L<5)
  off <- which(L <=1.5)
  pks <- which(L == pk)
  labeled <- c(pks,low,off)
  
  L[pks] <- "peak"
  L[low] <- "weak"
  L[off] <- "OFF"
  L[-labeled] <- "ON"
  COFG[i, 9:12] <- L
}
write.table(COFG, "/home/layal/OSTEO_2017/DESeq/On_OFF.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#141
COFG_pro <- COFG[which(COFG$Biotype == "protein_coding"),]
#96
COFG_noncod <- COFG[which(COFG$Biotype != "protein_coding"),]
#33
COFG_linc <- COFG[which(COFG$Biotype == "lincRNA"),]

#---------------------------------------------- separate biotype-----------------------
pro_type <- function(x){
  
  pro <- x[which(x$Gene.type == "protein_coding"),]
  return(pro)
}

non_type <- function(x){
  noncod <- x[which(x$Gene.type != "protein_coding"),]
  return(noncod)
}

linc_type <- function(x){ 
  linc <- x[which(x$Gene.type == "lincRNA"),]
  return(linc)
}
#---------------0-------------------
On0off3_pro <- pro_type(On0off3)
On0off3_non <- non_type(On0off3)
On0off3_linc <- linc_type(On0off3)

On0off6_pro <- pro_type(On0off6)
On0off6_non <- non_type(On0off6)
On0off6_linc <- linc_type(On0off6)

On0off12_pro <- pro_type(On0off12)
On0off12_non <- non_type(On0off12)
On0off12_linc <- linc_type(On0off12)
#----------------3------------------------
On3off0_pro <- pro_type(On3off0)
On3off0_non <- non_type(On3off0)
On3off0_linc <- linc_type(On3off0)

On3off6_pro <- pro_type(On3off6)
On3off6_non <- non_type(On3off6)
On3off6_linc <- linc_type(On3off6)

On3off12_pro <- pro_type(On3off12)
On3off12_non <- non_type(On3off12)
On3off12_linc <- linc_type(On3off12)

#----------------6----------------------------
On6off0_pro <- pro_type(On6off0)
On6off0_non <- non_type(On6off0)
On6off0_linc <- linc_type(On6off0)

On6off3_pro <- pro_type(On6off3)
On6off3_non <- non_type(On6off3)
On6off3_linc <- linc_type(On6off3)

On6off12_pro <- pro_type(On6off12)
On6off12_non <- non_type(On6off12)
On6off12_linc <- linc_type(On6off12)

#--------------------12 ------------------------
On12off0_pro <- pro_type(On12off0)
On12off0_non <- non_type(On12off0)
On12off0_linc <- linc_type(On12off0)

On12off3_pro <- pro_type(On12off3)
On12off3_non <- non_type(On12off3)
On12off3_linc <- linc_type(On12off3)

On12off6_pro <- pro_type(On12off6)
On12off6_non <- non_type(On12off6)
On12off6_linc <- linc_type(On12off6)

#==================================================plot counts===============================================================
runx <- plotCounts(ddseq, gene= "Runx2", intgroup = "condition", normalized = TRUE, transform = FALSE, returnData = TRUE)
runx$condition
runxN<- NormSF_Counts["Runx2", ]
runxMean <- c(mean(runxN[1:3]), mean(runxN[4:6]),mean(runxN[7:9]), mean(runxN[10:12]))
ggplot(runx, aes(x=condition, y=count)) +   geom_point(position=position_jitter(w=0.1,h=0)) 
