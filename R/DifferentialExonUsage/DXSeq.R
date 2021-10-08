# Author Layal
# version 1.1.0 last modified April20 2017
# Differential Exons usage between RNA Seq Osteoblast time points
# Input data: count files mapped with tophat to the reference GRCm38 release 86
# annotation file is from ensemble GRCm38 release 86



source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("DEXSeq")

library("DEXSeq")

#======================================================================== D12 - D0 =========================================================================================
#----------------count Exons ----------------------------------
# counting reads directory
inDirC = "/home/layal/OSTEO_2017/DEXseq/countsDXseq12_0"

#annotation directory
inDirG = "/home/layal/OSTEO_2017/DEXseq"

countFiles12_0 = list.files (inDirC, pattern=".txt$", full.names = TRUE)
basename(countFiles12_0)

flattenedFile = list.files (inDirG, pattern = "gff$", full.names = TRUE)
basename(flattenedFile)

sampleTable12_0 = data.frame(
  row.names = c("DayZero_R1","DayZero_R2","DayZero_R3" ,
                "DayTwelve_R1", "DayTwelve_R2", "DayTwelve_R3"),
  condition = c(rep("DayZero", 3),rep("DayTwelve",3)) )

dxd12_0 = DEXSeqDataSetFromHTSeq(
  countFiles12_0,
  sampleData=sampleTable12_0,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

colodxd12_0 <- colData(dxd12_0)

head( counts(dxd12_0), 5 )
#the count belonging to the exonic regions (’this’) 
head( featureCounts(dxd12_0), 5 )
dxd_countsRaw12_0 <- featureCounts(dxd12_0) 
dxf12_0 <- dxd12_0[ rowSums(featureCounts(dxd12_0)) > 1, ]

#details on the counting bins, list of the transcripts contain this exon and genomic coordinates
head( rowRanges(dxd12_0), 3 )
sampleAnnotation( dxd12_0) 

#--------------------Normalization ----------------------
#SMALLsubset

dxdn12_0 = estimateSizeFactors( dxf12_0 )

#-------------------dispersion estimation ----------------------
dxsd12_0 <- estimateDispersions( dxdn12_0 )

dxDEU12_0 = testForDEU( dxsd12_0 )

dxdFCZ = estimateExonFoldChanges( dxDEU12_0, fitExpToVar="condition", denominator = "DayZero")
dxr12_0 = DEXSeqResults( dxdFCZ )
DXR12_0 <-as.data.frame(dxr12_0)


table ( dxr12_0$padj < 0.01 )
table ( dxr12_0$padj < 0.001 )

mcols(dxr12_0)$description
# [1] "group/gene identifier"                                        "feature/exon identifier"                                     
# [3] "mean of the counts across samples in each feature/exon"       "exon dispersion estimate"                                    
# [5] "LRT statistic: full vs reduced"                               "LRT p-value: full vs reduced"                                
# [7] "BH adjusted p-values"                                         "exon usage coefficient"                                      
# [9] "exon usage coefficient"                                       "relative exon usage fold change"                             
# [11] "GRanges object of the coordinates of the exon/feature"        "matrix of integer counts, of each column containing a sample"
# [13] "list of transcripts overlapping with the exon"   

plotDEXSeq( dxr12_0, "ENSMUSG00000000194", expression=FALSE,displayTranscripts=TRUE,splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1, lwd=2 )

pdj <- 0.01
LFC <- log2(4)
sigExon <- DXR12_0[which(DXR12_0$padj <pdj & DXR12_0$exonBaseMean > 5), ]
DEU <- sigExon[ which(abs(sigExon$log2fold_DayTwelve_DayZero) > LFC) ,]

biomart <- read.table(paste(inDirG,"BioMart_GRCm38R86.txt",sep="/"), sep = "\t", header = TRUE)


library(hash)
gene_name_hash <- hash(keys= biomart$Ensembl.Gene.ID, values= biomart$Associated.Gene.Name)
biotype_hash <- hash(keys= biomart$Ensembl.Gene.ID, values= biomart$Gene.type)

DEUx <- DEU
for (i in 1:length(DEUx$groupID))
{
  if (grepl("\\+",DEUx$groupID[i]))
  {
    z<- unlist(strsplit(DEUx$groupID[i], "\\+"))
    x<- vector(mode = "character")
    y<- vector(mode = "character")

    for(k in 1:length(z))
    {
      if (is.null(gene_name_hash[[ z[k] ]])){
        x[k]<- "Not Available"
        y[k] <- "NA"
      }
      else{
        x[k] <- as.character(gene_name_hash[[ z[k] ]])
        y[k] <-as.character(biotype_hash[[ z[k] ]])
      }
    }
    DEUx$gene_name[i] <- paste(x ,sep = " + ", collapse = "+")
    DEUx$gene_t[i] <- paste(y ,sep = " + ", collapse = "+")
    
  }
  else
  {
    if (is.null(gene_name_hash[[ DEU$groupID[i] ]])){
      DEUx$gene_name[i] <- "Not Available"
      DEUx$gene_t[i] <- "NA"
    }
    else{
      DEUx$gene_name[i] <- as.character(gene_name_hash[[ DEUx$groupID[i] ]])
      DEUx$gene_t[i] <- as.character(biotype_hash[[ DEUx$groupID[i] ]])
    }
      
  }
}

write.table(DEUx, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DEU12_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.csv(DEUx ,"/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DEU12_0.csv",  row.names = TRUE)
# report
installed.packages("hwriter")
library("hwriter")
DEXSeqHTML( dxr, FDR=0.01 )

#======================================================================== D6 - D0 =========================================================================================
#----------------count Exons ----------------------------------
# counting reads directory
inDirC6_0 = "/home/layal/OSTEO_2017/DEXseq/countsDXseq6_0"



countFiles6_0 = list.files (inDirC6_0, pattern=".txt$", full.names = TRUE)
basename(countFiles6_0)

sampleTable6_0 = data.frame(
  row.names = c("DayZero_R1","DayZero_R2","DayZero_R3" ,
                "DaySix_R1", "DaySix_R2", "DaySix_R3"),
  condition = c(rep("DayZero", 3),rep("DaySix",3)) )

dxd6_0 = DEXSeqDataSetFromHTSeq(
  countFiles6_0,
  sampleData=sampleTable6_0,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

colData(dxd6_0)

head( counts(dxd6_0), 5 )
#the count belonging to the exonic regions (this) 
head( featureCounts(dxd6_0), 5 )
dxd_countsRaw6_0 <- featureCounts(dxd6_0) 
dxf6_0 <- dxd6_0[ rowSums(featureCounts(dxd6_0)) > 1, ]

#details on the counting bins, list of the transcripts contain this exon and genomic coordinates
head( rowRanges(dxd6_0), 3 )


#--------------------Normalization ----------------------
#SMALLsubset

dxdn6_0 = estimateSizeFactors( dxf6_0 )

#-------------------dispersion estimation ----------------------
dxsd6_0 <- estimateDispersions( dxdn6_0 )

dxDEU6 = testForDEU( dxsd6_0 )

dxdFC6_0 = estimateExonFoldChanges( dxDEU6, fitExpToVar="condition", denominator = "DayZero")
dxr6_0 = DEXSeqResults( dxdFC6_0 )
DXR6_0 <-as.data.frame(dxr6_0)

table ( dxr6_0$padj < 0.01 )
table ( dxr6_0$padj < 0.001 )

mcols(dxr6_0)$description

pdj <- 0.01
LFC <- log2(4)
sigExon6_0 <- DXR6_0[which(DXR6_0$padj <pdj & DXR6_0$exonBaseMean > 5), ]
DEU6_0 <- sigExon6_0[ which(abs(sigExon6_0$log2fold_DaySix_DayZero) > LFC) ,]

DEUx6_0 <- DEU6_0
for (i in 1:length(DEUx6_0$groupID))
{
  if (grepl("\\+",DEUx6_0$groupID[i]))
  {
    z<- unlist(strsplit(DEUx6_0$groupID[i], "\\+"))
    x<- vector(mode = "character")
    y<- vector(mode = "character")
    
    for(k in 1:length(z))
    {
      if (is.null(gene_name_hash[[ z[k] ]])){
        x[k]<- "Not Available"
        y[k] <- "NA"
      }
      else{
        x[k] <- as.character(gene_name_hash[[ z[k] ]])
        y[k] <-as.character(biotype_hash[[ z[k] ]])
      }
    }
      DEUx6_0$gene_name[i] <- paste(x ,sep = " + ", collapse = "+")
      DEUx6_0$gene_t[i] <-  paste(y ,sep = " + ", collapse = "+")
  }
  else
    {
      if (is.null(gene_name_hash[[ DEUx6_0$groupID[i] ]])){
        DEUx6_0$gene_name[i] <- "Not Available"
        DEUx6_0$gene_t[i] <- "NA"
      }
      else{
        DEUx6_0$gene_name[i] <- as.character(gene_name_hash[[ DEUx6_0$groupID[i] ]])
        DEUx6_0$gene_t[i] <- as.character(biotype_hash[[ DEUx6_0$groupID[i] ]])
        
      }
        
    }
}

write.table(DEUx6_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DEU6_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
#======================================================================== D3 - D0 =========================================================================================
#----------------count Exons ----------------------------------
# counting reads directory
inDirC3_0 = "/home/layal/OSTEO_2017/DEXseq/countsDXseq3_0"



countFiles3_0 = list.files (inDirC3_0, pattern=".txt$", full.names = TRUE)
basename(countFiles3_0)

sampleTable3_0 = data.frame(
  row.names = c("DayZero_R1","DayZero_R2","DayZero_R3" ,
                "DayThree_R1", "DayThree_R2", "DayThree_R3"),
  condition = c(rep("DayZero", 3),rep("DayThree",3)) )

dxd3_0 = DEXSeqDataSetFromHTSeq(
  countFiles3_0,
  sampleData=sampleTable3_0,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

colData(dxd3_0)

head( counts(dxd3_0), 5 )
#the count belonging to the exonic regions (this) 
head( featureCounts(dxd3_0), 5 )
dxd_countsRaw3_0 <- featureCounts(dxd3_0) 
dxf3_0 <- dxd3_0[ rowSums(featureCounts(dxd3_0)) > 1, ]

#details on the counting bins, list of the transcripts contain this exon and genomic coordinates
head( rowRanges(dxd3_0), 3 )


#--------------------Normalization ----------------------
#SMALLsubset

dxdn3_0 = estimateSizeFactors( dxf3_0 )

#-------------------dispersion estimation ----------------------
dxsd3_0 <- estimateDispersions( dxdn3_0 )

dxDEU3_0 = testForDEU( dxsd3_0 )

dxdFC3_0 = estimateExonFoldChanges( dxDEU3_0, fitExpToVar="condition", denominator = "DayZero")
dxr3_0 = DEXSeqResults( dxdFC3_0 )
DXR3_0 <-as.data.frame(dxr3_0)

table ( dxr3_0$padj < 0.01 )
table ( dxr3_0$padj < 0.001 )

mcols(dxr3_0)$description

pdj <- 0.01
LFC <- log2(4)
sigExon3_0 <- DXR3_0[which(DXR3_0$padj <pdj & DXR3_0$exonBaseMean > 5), ]
DEU3_0 <- sigExon3_0[ which(abs(sigExon3_0$log2fold_DayThree_DayZero) > LFC) ,]

DEUx3_0 <- DEU3_0

for (i in 1:length(DEUx3_0$groupID))
{
  if (grepl("\\+",DEUx3_0$groupID[i]))
  {
    z<- unlist(strsplit(DEUx3_0$groupID[i], "\\+"))
    x<- vector(mode = "character")
    y<- vector(mode = "character")
    
    for(k in 1:length(z))
    {
      if (is.null(gene_name_hash[[ z[k] ]])){
        x[k]<- "Not Available"
        y[k] <- "NA"
      }
      else{
        x[k] <- as.character(gene_name_hash[[ z[k] ]])
        y[k] <-as.character(biotype_hash[[ z[k] ]])
      }
    }
      DEUx3_0$gene_name[i] <- paste(x ,sep = " + ", collapse = "+")
      DEUx3_0$gene_t[i] <-  paste(y ,sep = " + ", collapse = "+")
  }
  else
    {
      if (is.null(gene_name_hash[[ DEUx3_0$groupID[i] ]])){
        DEUx3_0$gene_name[i] <- "Not Available"
        DEUx3_0$gene_t[i] <- "NA"
      }
      else{
        DEUx3_0$gene_name[i] <- as.character(gene_name_hash[[ DEUx3_0$groupID[i] ]])
        DEUx3_0$gene_t[i] <- as.character(biotype_hash[[ DEUx3_0$groupID[i] ]])
        
      }
      
    }
}

write.table(DEUx3_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DEU3_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#======================================================================== D12 - D3 =========================================================================================
#----------------count Exons ----------------------------------
# counting reads directory
inDirC12_3 = "/home/layal/OSTEO_2017/DEXseq/countsDXseq12_3"



countFiles12_3 = list.files (inDirC12_3, pattern=".txt$", full.names = TRUE)
basename(countFiles12_3)

sampleTable12_3 = data.frame(
  row.names = c("DayThree_R1", "DayThree_R2", "DayThree_R3",
                "DayTwelve_R1", "DayTwelve_R2", "DayTwelve_R3"),
  condition = c(rep("DayThree",3), rep("DayTwelve", 3)))

dxd12_3 = DEXSeqDataSetFromHTSeq(
  countFiles12_3,
  sampleData=sampleTable12_3,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

colData(dxd12_3)

head( counts(dxd12_3), 5 )
#the count belonging to the exonic regions (’this’) 
head( featureCounts(dxd12_3), 5 )
dxd_countsRaw12_3 <- featureCounts(dxd12_3) 
dxf12_3 <- dxd12_3[ rowSums(featureCounts(dxd12_3)) > 1, ]

#details on the counting bins, list of the transcripts contain this exon and genomic coordinates
head( rowRanges(dxd12_3), 3 )


#--------------------Normalization ----------------------
#SMALLsubset

dxdn12_3 = estimateSizeFactors( dxf12_3 )

#-------------------dispersion estimation ----------------------
dxsd12_3 <- estimateDispersions( dxdn12_3 )

dxDEU12_3 = testForDEU( dxsd12_3 )

dxdFC12_3 = estimateExonFoldChanges( dxDEU12_3, fitExpToVar="condition", denominator = "DayThree")
dxr12_3 = DEXSeqResults( dxdFC12_3 )
DXR12_3 <-as.data.frame(dxr12_3)

table ( dxr12_3$padj < 0.01 )
table ( dxr12_3$padj < 0.001 )

mcols(dxr12_3)$description

pdj <- 0.01
LFC <- log2(3)
sigExon12_3 <- DXR12_3[which(DXR12_3$padj <pdj & DXR12_3$exonBaseMean > 5), ]
DEU12_3 <- sigExon12_3[ which(abs(sigExon12_3$log2fold_DayTwelve_DayThree) > LFC) ,]

DEUx12_3 <- DEU12_3

for (i in 1:length(DEUx12_3$groupID))
{
  if (grepl("\\+",DEUx12_3$groupID[i]))
  {
    z<- unlist(strsplit(DEUx12_3$groupID[i], "\\+"))
    x<- vector(mode = "character")
    y<- vector(mode = "character")
    
    for(k in 1:length(z))
    {
      if (is.null(gene_name_hash[[ z[k] ]])){
        x[k]<- "Not Available"
        y[k] <- "NA"
      }
      else{
        x[k] <- as.character(gene_name_hash[[ z[k] ]])
        y[k] <- as.character(biotype_hash[[ z[k] ]])
      }
    }
      DEUx12_3$gene_name[i] <- paste(x ,sep = " + ", collapse = "+")
      DEUx12_3$gene_t[i] <-  paste(y ,sep = " + ", collapse = "+")
      
  }
  else
    {
      if (is.null(gene_name_hash[[ DEUx12_3$groupID[i] ]])){
        DEUx12_3$gene_name[i] <- "Not Available"
        DEUx12_3$gene_t[i] <- "NA"
      }
      else{
        DEUx12_3$gene_name[i] <- as.character(gene_name_hash[[ DEUx12_3$groupID[i] ]])
        DEUx12_3$gene_t[i] <- as.character(biotype_hash[[ DEUx12_3$groupID[i] ]])
        
      }
      
    }
}

write.table(DEUx12_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DEU12_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
#======================================================================== D12 - D6 =========================================================================================
#----------------count Exons ----------------------------------
# counting reads directory
inDirC12_6 = "/home/layal/OSTEO_2017/DEXseq/countsDXseq12_6"



countFiles12_6 = list.files (inDirC12_6, pattern=".txt$", full.names = TRUE)
basename(countFiles12_6)

sampleTable12_6 = data.frame(
  row.names = c("DaySix_R1", "DaySix_R2", "DaySix_R3",
                "DayTwelve_R1", "DayTwelve_R2", "DayTwelve_R3"),
  condition = c(rep("DaySix",3), rep("DayTwelve", 3)))

dxd12_6 = DEXSeqDataSetFromHTSeq(
  countFiles12_6,
  sampleData=sampleTable12_6,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

colData(dxd12_6)

head( counts(dxd12_6), 5 )
#the count belonging to the exonic regions (’this’) 
head( featureCounts(dxd12_6), 5 )
dxd_countsRaw12_6 <- featureCounts(dxd12_6) 
dxf12_6 <- dxd12_6[ rowSums(featureCounts(dxd12_6)) > 1, ]

#details on the counting bins, list of the transcripts contain this exon and genomic coordinates
head( rowRanges(dxd12_6), 3 )


#--------------------Normalization ----------------------
#SMALLsubset

dxdn12_6 = estimateSizeFactors( dxf12_6 )

#-------------------dispersion estimation ----------------------
dxsd12_6 <- estimateDispersions( dxdn12_6 )

dxDEU12_6 = testForDEU( dxsd12_6 )

dxdFC12_6 = estimateExonFoldChanges( dxDEU12_6, fitExpToVar="condition", denominator = "DaySix")
dxr12_6 = DEXSeqResults( dxdFC12_6 )
DXR12_6 <-as.data.frame(dxr12_6)

table ( dxr12_6$padj < 0.01 )
table ( dxr12_6$padj < 0.001 )

mcols(dxr12_6)$description

pdj <- 0.01
LFC <- log2(3)

sigExon12_6 <- DXR12_6[which(DXR12_6$padj <pdj & DXR12_6$exonBaseMean > 5), ]
DEU12_6 <- sigExon12_6[ which(abs(sigExon12_6$log2fold_DayTwelve_DaySix) > LFC) ,]

DEUx12_6 <- DEU12_6

for (i in 1:length(DEUx12_6$groupID))
{
  if (grepl("\\+",DEUx12_6$groupID[i]))
  {
    z<- unlist(strsplit(DEUx12_6$groupID[i], "\\+"))
    x<- vector(mode = "character")
    y<- vector(mode = "character")
    
    for(k in 1:length(z))
    {
      if (is.null(gene_name_hash[[ z[k] ]])){
        x[k]<- "Not Available"
        y[k] <- "NA"
      }
      else{
        x[k] <-as.character(gene_name_hash[[ z[k] ]])
        y[k] <- as.character(biotype_hash[[ z[k] ]])
      }
    }
      DEUx12_6$gene_name[i] <- paste(x ,sep = " + ", collapse = "+")
      DEUx12_6$gene_t[i] <-  paste(y ,sep = " + ", collapse = "+")
  }
  else
    {
      if (is.null(gene_name_hash[[ DEUx12_6$groupID[i] ]])){
        DEUx12_6$gene_name[i] <- "Not Available"
        DEUx12_6$gene_t[i] <- "NA"
        }
      else{
        DEUx12_6$gene_name[i] <- as.character(gene_name_hash[[ DEUx12_6$groupID[i] ]])
        DEUx12_6$gene_t[i] <- as.character(biotype_hash[[ DEUx12_6$groupID[i] ]])
        }
      }
}

write.table(DEUx12_6, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DEU12_6.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#======================================================================== D6 - D3 =========================================================================================
#----------------count Exons ----------------------------------
# counting reads directory
inDirC6_3 = "/home/layal/OSTEO_2017/DEXseq/countsDXseq6_3"



countFiles6_3 = list.files (inDirC6_3, pattern=".txt$", full.names = TRUE)
basename(countFiles6_3)

sampleTable6_3 = data.frame(
  row.names = c("DayThree_R1", "DayThree_R2", "DayThree_R3","DaySix_R1", "DaySix_R2", "DaySix_R3"),
  condition = c(rep("DayThree",3), rep("DaySix", 3)))

dxd6_3 = DEXSeqDataSetFromHTSeq(
  countFiles6_3,
  sampleData=sampleTable6_3,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

colData(dxd6_3)

head( counts(dxd6_3), 5 )
#the count belonging to the exonic regions ("this") 
head( featureCounts(dxd6_3), 5 )
dxd_countsRaw6_3 <- featureCounts(dxd6_3) 
dxf6_3 <- dxd6_3[ rowSums(featureCounts(dxd6_3)) > 1, ]

#details on the counting bins, list of the transcripts contain this exon and genomic coordinates
head( rowRanges(dxd6_3), 3 )


#--------------------Normalization ----------------------
#SMALLsubset

dxdn6_3 = estimateSizeFactors( dxf6_3 )

#-------------------dispersion estimation ----------------------
dxsd6_3 <- estimateDispersions( dxdn6_3 )

dxDEU6_3 = testForDEU( dxsd6_3 )

dxdFC6_3 = estimateExonFoldChanges( dxDEU6_3, fitExpToVar="condition", denominator = "DayThree")
dxr6_3 = DEXSeqResults( dxdFC6_3 )
DXR6_3 <-as.data.frame(dxr6_3)

table ( dxr6_3$padj < 0.01 )
table ( dxr6_3$padj < 0.001 )

mcols(dxr6_3)$description

pdj <- 0.01
LFC <- log2(3)

sigExon6_3 <- DXR6_3[which(DXR6_3$padj <pdj & DXR6_3$exonBaseMean > 5), ]
DEU6_3 <- sigExon6_3[ which(abs(sigExon6_3$log2fold_DaySix_DayThree) > LFC) ,]

DEUx6_3 <- DEU6_3

for (i in 1:length(DEUx6_3$groupID))
{
  if (grepl("\\+",DEUx6_3$groupID[i]))
  {
    z<- unlist(strsplit(DEUx6_3$groupID[i], "\\+"))
    x<- vector(mode = "character")
    y<- vector(mode = "character")
    
    for(k in 1:length(z))
    {
      if (is.null(gene_name_hash[[ z[k] ]])){
        x[k]<- "Not Available"
        y[k] <- "NA"
      }
      else{
        x[k] <-as.character(gene_name_hash[[ z[k] ]])
        y[k] <- as.character(biotype_hash[[ z[k] ]])
      }
    }
      DEUx6_3$gene_name[i] <- paste(x ,sep = " + ", collapse = "+")
      DEUx6_3$gene_t[i] <-  paste(y ,sep = " + ", collapse = "+")
      
  }
  else
  {
    if (is.null(gene_name_hash[[ DEUx12_6$groupID[i] ]])){
      DEUx6_3$gene_name[i] <- "Not Available"
      DEUx6_3$gene_t[i] <- "NA"
    }
    else{
      DEUx6_3$gene_name[i] <- as.character(gene_name_hash[[ DEUx12_6$groupID[i] ]])
      DEUx6_3$gene_t[i] <- as.character(biotype_hash[[ DEUx12_6$groupID[i] ]])
      
    }
    
  }
}

write.table(DEUx6_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DEU6_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#=============================================================================================================================================================================
#=========================================== Choose only the protein transcript==================================================
# this function get the name and the type of the transcripts and save them as a separate vector

# we passed DX is the DXseq object object , and TX is a vectore that will store all the transcripts

# read biomart table with transcript type

trans_biomart <- read.table(paste(inDirG,"transcriptsBiotype.txt",sep="/"), sep = "\t", header = TRUE)

transName_hash <- hash(keys= trans_biomart$Ensembl.Transcript.ID, values= trans_biomart$Associated.Transcript.Name)
transBiotype_hash <- hash(keys= trans_biomart$Ensembl.Transcript.ID, values= trans_biomart$Transcript.type)

trans_ProCod <- function(DX, TX){
  t <- as.character(DX$transcripts)
  
  for (i in 1:length(t))
  {
    if (grepl("c",t[i]))
    {
      c <-gsub( 'c','', (gsub("[[:punct:]]",'', t[i])))
      z<- unlist(strsplit(c, " "))
      x<- vector(mode = "character")
      y<- vector(mode = "character")
      for(k in 1:length(z))
      {
        if (is.null(transName_hash[[ z[k] ]])){
          x[k]<- "Not Available"
          y[k] <- "NA"
        }
        else{
          x[k] <-as.character(transName_hash[[ z[k] ]])
          y[k] <- as.character(transBiotype_hash[[ z[k] ]])
        }
       
      }
      DX$transName[i] <- paste(x ,sep = "  +  ", collapse = "+")
      DX$transBiotype[i] <-  paste(y ,sep = "  +  ", collapse = "+")
      TX <- append(TX, z)
    }
    else{
      TX <- append(TX, t[i])
      if (is.null(transName_hash[[ t[i] ]])){
        DX$transName[i] <- "Not Available"
        DX$transBiotype[i] <- "NA"
      }
      else{
        DX$transName[i] <- as.character(transName_hash[[ t[i] ]])
        DX$transBiotype[i] <- as.character(transBiotype_hash[[ t[i] ]])
      }
    }
  }
  
  
  DXPro <- DX[ which(grepl("protein_coding", DX$transBiotype)) ,c(1,23,24,22,25,26,2,3,7,8,9,10,12,13,14,15)]
  
  resPro = list(DX = DXPro, Trans= TX)
  return(resPro)
}

Trans12_0 <- vector(mode = "character")
DXPro12_0 <- trans_ProCod(DEUx, Trans12_0)$DX
Trans12_0 <- trans_ProCod(DEUx, Trans12_0)$Trans

Trans12_3 <- vector(mode = "character")
DXPro12_3 <- trans_ProCod(DEUx12_3, Trans12_3)$DX
Trans12_3 <- trans_ProCod(DEUx12_3, Trans12_3)$Trans

Trans12_6 <- vector(mode = "character")
DXPro12_6 <- trans_ProCod(DEUx12_6, Trans12_6)$DX
Trans12_6 <- trans_ProCod(DEUx12_6, Trans12_6)$Trans

Trans6_0 <- vector(mode = "character")
DXPro6_0 <- trans_ProCod(DEUx6_0, Trans6_0)$DX
Trans6_0 <- trans_ProCod(DEUx6_0, Trans6_0)$Trans


Trans6_3 <- vector(mode = "character")
DXPro6_3 <- trans_ProCod(DEUx6_3, Trans6_3)$DX
Trans6_3 <- trans_ProCod(DEUx6_3, Trans6_3)$Trans

Trans3_0 <- vector(mode = "character")
DXPro3_0 <- trans_ProCod(DEUx3_0, Trans3_0)$DX
Trans3_0 <- trans_ProCod(DEUx3_0, Trans3_0)$Trans

write.table(DXPro12_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DXPro12_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DXPro12_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DXPro12_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DXPro12_6, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DXPro12_6.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DXPro6_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DXPro6_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DXPro3_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DXPro3_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DXPro6_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DXPro6_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#==========================================================================================================================================================
#================================================ function for exon coding========================================================================

# in this function we need to see if the exon bin (in DXSeq object), it is not the full exon, most of the time is part of exon,
# so from exons features from BioMart in Ensemble, 
#    Exon.Chr.Start..bp.   Exon.Chr.End..bp.   Ensembl.Transcript.ID   Genomic.coding.start    Genomic.coding.end    Strand
#        15356               15422             ENSMUST00000082423            NA                      NA                -1
#        14145               15288            ENSMUST00000082421            14145                   15288               1


# our parameters
# S : is the start of the full exon (Exon.Chr.Start..bp.)
# E : the end of the full exon 
# Sc : Genomic.coding.start
# Ec : Genomic.coding.end
# Sb : is the exon bin start, from Dxseq table genomicData.start
# Eb : is the end of exon bin

ExonsCor <-  read.table(paste(inDirG,"exonsCor.txt",sep="/"), sep = "\t", header = TRUE)

Trans <- Reduce(union, list(Trans12_0,Trans12_3,Trans12_6, Trans3_0,Trans6_0,Trans6_3))
exonCor <- ExonsCor[which(ExonsCor$Ensembl.Transcript.ID %in% Trans),]
CodStart_hash <- hash(keys= exonCor$Exon.Chr.Start..bp. , values=exonCor$Genomic.coding.start)
CodEnd_hash <- hash(keys= exonCor$Exon.Chr.End..bp., values=exonCor$Genomic.coding.end)


ExonORF <- function(DX, trans){
  BMexon <- ExonsCor[which(ExonsCor$Ensembl.Transcript.ID %in% trans),]
  for(j in 1:nrow(DX)){
    Sb <- DX$genomicData.start[j]
    
    Eb <- DX$genomicData.end[j]
   
    
    for(i in 1:nrow(BMexon)){
      S <- BMexon$Exon.Chr.Start..bp.[i]
      E <- BMexon$Exon.Chr.End..bp.[i]
      Sc <- BMexon$Genomic.coding.start[i]
      Ec <- BMexon$Genomic.coding.end[i]
      
      if(Sb>=S && Eb<=E){
        cat("\n start of bin Sb and the end is inside S,E ")
        
        
        if(!is.na(Sc) && !is.na(Ec)){
          cat("\n start coding is integer ")
          
          if(Sb>=Sc && Eb<=Ec){
            cat("\n  exonic bin is coding  ")
            DX$ExonType[j] <- "coding_Exon"
          }
          else{
            print("********************NOT ORF ")
            DX$ExonType[j] <- "out_ORF"
          }
        }
        else{
          print(" ****************** Sc is NA : exon NOT coding ")
          DX$ExonType[j] <- "NOT_Coding"
        }
      }
     
      
    }
  }
  
  return(DX)
}

DX12_0 <- ExonORF(DXPro12_0, Trans12_0)
DX12_3 <- ExonORF(DXPro12_3, Trans12_3)
DX12_6 <- ExonORF(DXPro12_6, Trans12_6)
DX6_0 <- ExonORF(DXPro6_0, Trans6_0)
DX3_0 <- ExonORF(DXPro3_0, Trans3_0)
DX6_3 <- ExonORF(DXPro6_3, Trans6_3)
write.table(DX12_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DX12_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DX12_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DX12_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DX12_6, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DX12_6.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DX6_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DX6_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DX3_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DX3_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DX6_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/DX6_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

ExCod12_0 <- DX12_0[which(DX12_0$ExonType == "coding_Exon"), ]
ExCod12_3 <- DX12_3[which(DX12_3$ExonType == "coding_Exon"), ]
ExCod12_6 <- DX12_6[which(DX12_6$ExonType == "coding_Exon"), ]
ExCod6_0 <- DX6_0[which(DX6_0$ExonType == "coding_Exon"), ]
ExCod3_0 <- DX3_0[which(DX3_0$ExonType == "coding_Exon"), ]
ExCod6_3 <- DX6_3[which(DX6_3$ExonType == "coding_Exon"), ]

write.table(ExCod12_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/ExCod12_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ExCod12_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/ExCod12_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ExCod12_6, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/ExCod12_6.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ExCod6_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/ExCod6_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ExCod3_0, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/ExCod3_0.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ExCod6_3, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DXSeq/output_July/ExCod6_3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)



