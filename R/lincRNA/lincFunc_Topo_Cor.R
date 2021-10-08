#This script to get the related RNAlinc to gene coding
# Author Layal v.0.0.2 in 7June 2017

library(DESeq2)
library(hash)
library(gplots)

pathIn="/home/layal/OSTEO_2017/DESeq"
path="/home/layal/OSTEO_2017/lnc_RNA"

#=============================================== creat DESeq object =====================================================
# get the counts matrix
FC_countreads <- readRDS(paste(pathIn,"FC_counts.rds",sep = "/"))
countdata <- FC_countreads$counts
# change names
colnames(countdata) <- gsub(pattern = "X.home.layal.Solexa.layal.OSTEO.Align_Tophat_GRCm38.|.accepted_hits.bam",replacement = "",x = colnames(countdata))
New_names <- c("DayZero_R1","DayZero_R2","DayZero_R3" ,"DayThree_R1","DayThree_R2","DayThree_R3","DaySix_R1","DaySix_R2","DaySix_R3", "DayTwelve_R1", "DayTwelve_R2", "DayTwelve_R3")
colnames(countdata) <- New_names

# creat DESeqDataSet Object
condition <- factor(c(rep("DayZero", 3),rep("DayThree",3),rep("DaySix",3),rep("DayTwelve",3)))
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
# Remove the counts less than 1
ddsf <- dds[ rowSums(counts(dds)) > 1, ]
# get the DESeq object
ddseq <- DESeq(ddsf)

#=========================================================== Results ======================================================

# Day 12 to Day 0 comparison :

resTwelve_Zero <- results(ddseq, contrast=c("condition", "DayTwelve", "DayZero"),pAdjustMethod = "bonferroni", alpha = 0.01,
                          cooksCutoff = TRUE, independentFiltering = TRUE)
res12_0 <- as.data.frame(resTwelve_Zero)

resSix_Zero <- results(ddseq, contrast=c("condition", "DaySix" ,"DayZero") , pAdjustMethod = "bonferroni", alpha = 0.01,
                       cooksCutoff = TRUE, independentFiltering = TRUE)
res6_0 <- as.data.frame(resSix_Zero)

resThree_Zero <- results(ddseq, contrast=c("condition", "DayThree", "DayZero"), pAdjustMethod = "bonferroni", alpha = 0.01,
                         cooksCutoff = TRUE, independentFiltering = TRUE)
res3_0 <- as.data.frame(resThree_Zero)

resTwelve_Three <- results(ddseq, contrast=c("condition", "DayTwelve", "DayThree"),pAdjustMethod = "bonferroni", alpha = 0.01,
                           cooksCutoff = TRUE, independentFiltering = TRUE)
res12_3 <- as.data.frame(resTwelve_Three)

resTwelve_Six <- results(ddseq, contrast=c("condition", "DayTwelve", "DaySix"),pAdjustMethod = "bonferroni", alpha = 0.01,
                         cooksCutoff = TRUE, independentFiltering = TRUE)
res12_6 <- as.data.frame(resTwelve_Six)

resSix_Three <- results(ddseq, contrast=c("condition", "DaySix", "DayThree"),pAdjustMethod = "bonferroni", alpha = 0.01,
                        cooksCutoff = TRUE, independentFiltering = TRUE)
res6_3 <- as.data.frame(resSix_Three)

# =============================== data frame with all comparisons ======================================================
# we need on data frame contain all comparison with logarithm fold change and adj pvalue
DifG <- data.frame(row.names = rownames(ddseq))

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
}

#======================================Significant Diffrential expressed genes ===================================================
# We need to get the significant differential expressed genes , minimum the pvalue 0.01 and at least there is 1.5 fold change  

pdj <-0.01
FC <- 1.5
#12932
SigDifGen <- DifGA[ which( (DifGA$Padj3_0 < pdj ) | (DifGA$Padj6_0 < pdj ) | (DifGA$Padj12_0 < pdj)
                           | (DifGA$Padj12_3 < pdj ) | (DifGA$Padj12_6 < pdj) | (DifGA$Padj6_3 < pdj)), ]
#11928
DifGE <- SigDifGen[which( (abs(SigDifGen$LFC3_0) > log2(FC)) | (abs(SigDifGen$LFC6_0)>log2(FC)) | (abs(SigDifGen$LFC12_0)>log2(FC))
                          | (abs(SigDifGen$LFC12_3)> log2(FC)) | (abs(SigDifGen$LFC12_6)>log2(FC)) | (abs(SigDifGen$LFC6_3) >log2(FC))),]   

#=================================== get the biotype and the coordinates ==============================================
# I installed the BioMart table from Ensembl archive R.86 contains the genomic cordinates, chromosome name and strand, and gene type..
BioType <- read.table(paste(path,"biomart86TopFun.txt", sep = "/") ,sep = "\t", header = T)


Diff_BioType <- subset(BioType, BioType$Associated.Gene.Name %in% (rownames(DifGE)))

# there are repeated genes because the number of genes in DifGE= 11928 while in Diff_BioType = 12009.
#so we need to delete the dublicates but need to save the genes which has two biotype as proteincoding and linc RNA

pureBioType <- as.data.frame( Diff_BioType[-which(duplicated(Diff_BioType$Associated.Gene.Name) == TRUE),])


Dup <- Diff_BioType[which(duplicated(Diff_BioType$Associated.Gene.Name) == TRUE),]
DupA <- subset(Diff_BioType, Diff_BioType$Associated.Gene.Name %in% Dup$Associated.Gene.Name)
tlinc<-DupA[which(DupA$Gene.type == "lincRNA"),]
tpro<-DupA[which(DupA$Gene.type == "protein_coding"),]

#  The single gene name in tpro(the protein cod from duplicated gene names DupA) .. this means the duplicated name has different biotype
td <-tpro[which(duplicated(tpro$Associated.Gene.Name) == TRUE),]
ts <- subset(tpro, tpro$Associated.Gene.Name %in% td$Associated.Gene.Name)

# the single genes without pair in protein coding
tm <- tpro[! tpro$Associated.Gene.Name %in% ts$Associated.Gene.Name, ]

#the gene we must check that there biotype is protein coding 
genestomod <- as.character(tm$Associated.Gene.Name)
pureBioType1 <- pureBioType

for (i in genestomod) {
  pureBioType1[which(pureBioType1$Associated.Gene.Name == i ),3]<-"protein_coding"
}

p1 <- pureBioType[which(pureBioType$Gene.type == "protein_coding"),]
p2 <- pureBioType1[which(pureBioType1$Gene.type == "protein_coding"),]


Biotype_DGE <- pureBioType1
#================================================= Building hash of biotype ========================================================================
hash_biotype <- hash(keys= Biotype_DGE$Associated.Gene.Name , values= Biotype_DGE$Gene.type)
hash_G.start <- hash(keys= Biotype_DGE$Associated.Gene.Name , values= Biotype_DGE$Gene.Start..bp.)
hash_G.end   <- hash(keys= Biotype_DGE$Associated.Gene.Name , values= Biotype_DGE$Gene.End..bp.)
hash_chr <- hash(keys= Biotype_DGE$Associated.Gene.Name, values=Biotype_DGE$Chromosome.Name)
hash_strand <- hash(keys= Biotype_DGE$Associated.Gene.Name, values= Biotype_DGE$Strand)
hash_Ensembl <- hash(keys= Biotype_DGE$Associated.Gene.Name , values= Biotype_DGE$Ensembl.Gene.ID)

# add the biotype and the genomics cordinates to the diffrential expression genes matrix
for (i in 1:length(rownames(DifGE))) {
  DifGE$GeneType[i] <- as.character( hash_biotype[[ rownames(DifGE)[i] ]])
  DifGE$Chromosome[i] <- as.character( hash_chr[[ rownames(DifGE)[i] ]])
  DifGE$Start[i] <- hash_G.start[[ rownames(DifGE)[i] ]]
  DifGE$End[i] <- hash_G.end[[ rownames(DifGE)[i] ]]
  DifGE$Strand[i] <- hash_strand[[ rownames(DifGE)[i] ]]
  DifGE$Ensembl[i] <- as.character(hash_Ensembl[[ rownames(DifGE)[i] ]])
  
}

#10778
DifGE_pro <- DifGE[which(DifGE$GeneType =="protein_coding"),]
#285
DifGE_linc <- DifGE[which(DifGE$GeneType =="lincRNA"),]


# ============================== topological relation and co-expression correlation algorithm =====================================================================
# we need the fold change values to calculate the correlation of each protein coding genes with link RNA.
# the correlation is calculated for the colums so we need to transpose the matris and take just the LFC od the comparisons
cor_DifPro <- t(DifGE_pro[,1:6])
cor_DifLinc <- t(DifGE_linc[,1:6])


# here the algorithm is working as the following:
# 1. we take the significat Diffirentially expressed genes DifGE. the protein coding "DifGE_pro" and the linRNA.
# 2. if the pro_cod gene and the lincRNA on the same chromosome, see in they are on the same strand, then check if they are in the neiborhood
#    we take the distance from lincRNA as variable D
# 3. if the protein coding gene has correlation > 0.9 and the pvalue of correlation test is < 0.01 , then consider the lincRNA related to this protein coding gene

D <- 2E6
cr <- 0.95
pv<- 0.01

#1407
Pro_Linc <- data.frame()
for(i in 1:length(rownames(DifGE_linc))){
  for(j in 1:length(rownames(DifGE_pro))){
    
    chr <- DifGE_linc$Chromosome[i]
    
    if ( chr== DifGE_pro$Chromosome[j] ){
      cat("\n  chr of lnc-RNA the same as chr of pro-cod\n")
      linc_Strand <-DifGE_linc$Strand[i]
      pro_Strand <- DifGE_pro$Strand[j]
      Start_Pro <-DifGE_pro$Start[j]
      End_Pro <- DifGE_pro$End[j]
      Start_Linc <-DifGE_linc$Start[i]
      End_Linc <-DifGE_linc$End[i]
      Protein_cod <- rownames(DifGE_pro[j,])
      LincRNA <-rownames(DifGE_linc[i,])  
     
       if( (Start_Pro %in% seq(Start_Linc - D , Start_Linc + D)) || (End_Pro %in% seq(End_Linc - D, End_Linc + D)) ){
        cat("\n\n ...............  pro-cod inside the range...........................\n\n")
         Correlation <- cor(cor_DifLinc[,i],cor_DifPro[,j])
         Pvalue <- cor.test(cor_DifLinc[,i],cor_DifPro[,j])$p.value
         if (abs(Correlation) > cr && Pvalue < pv){
           cat("\n\n .......................    correlation crateria achieved............................ \n\n")
           Pro_Linc = rbind(Pro_Linc, data.frame(LincRNA, Protein_cod, Correlation, Pvalue,chr, Start_Linc,End_Linc, linc_Strand,Start_Pro, End_Pro, pro_Strand))
          
         }
        
      }
      
    }
    
  }
}

write.table(Pro_Linc, "/home/layal/OSTEO_2017/lnc_RNA/PorCod_linc_allpairs.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

saveRDS(Pro_Linc,"/home/layal/OSTEO_2017/lnc_RNA/pro_linc0.95.rds")
#235
pure_linc <- Pro_Linc[-which(duplicated(Pro_Linc$LincRNA)==TRUE),]
#1252
pure_pro <- Pro_Linc[-which(duplicated(Pro_Linc$Protein_cod)==TRUE),]


#50
notLinc<- DifGE_linc[! rownames(DifGE_linc) %in% pure_linc$LincRNA,]


Lincs <- unique(Pro_Linc$LincRNA)

posCor <- data.frame()
negCor <- data.frame()
for(i in Lincs){
  prots <- Pro_Linc[which(Pro_Linc$LincRNA == i), ]
  if(max(prots$Correlation) > 0){
    posCor <- rbind(posCor, prots[which.max(prots$Correlation),])
    }
  if(min(prots$Correlation) < 0){
    negCor <- rbind(negCor, prots[which.min(prots$Correlation),])
    }
  
}

write.table(posCor, "/home/layal/OSTEO_2017/lnc_RNA/PorCod_linc_pos.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(negCor, "/home/layal/OSTEO_2017/lnc_RNA/PorCod_linc_neg.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#43
alwPosCoe <- posCor[ ! posCor$LincRNA %in% negCor$LincRNA,]

#28
alwNegCor <- negCor[ ! negCor$LincRNA %in% posCor$LincRNA,]
dpro <- posCor[ which(posCor$Protein_cod %in% negCor$Protein_cod),]
dnpro <- negCor[ which(negCor$Protein_cod %in% dpro$Protein_cod ),]

#============================================== Data Frame contain all the linc correlated ===================================
proPos_hash <-hash(keys= as.character(posCor$LincRNA) , values= as.character(posCor$Protein_cod))
St_proPos_hash <- hash(keys= as.character(posCor$LincRNA) , values= posCor$Start_Pro)
End_proPos_hash <- hash(keys= as.character(posCor$LincRNA) , values= posCor$End_Pro)

proNeg_hash <- hash(keys= as.character(negCor$LincRNA) , values= as.character(negCor$Protein_cod))
St_proNeg_hash <- hash(keys= as.character(negCor$LincRNA) , values= negCor$Start_Pro)
End_proNeg_hash <- hash(keys= as.character(negCor$LincRNA) , values= negCor$End_Pro)

linc_Pro_pairs <- data.frame(row.names = rownames(pure_linc))

for (i in 1:length(pure_linc$LincRNA)){
  Linc_RNA <- as.character(pure_linc$LincRNA[i])
  Chrom <- as.character (pure_linc$chr[i])
  Start_linc <- pure_linc$Start_Linc[i]
  End_linc <- pure_linc$End_Linc[i]
  Start_posPro <- St_proPos_hash [[Linc_RNA]]
  End_posPro <- End_proPos_hash[[Linc_RNA]]
  Start_negPro <- St_proNeg_hash [[Linc_RNA]]
  End_negPro <- End_proNeg_hash[[Linc_RNA]]
       
  if(is.null(proPos_hash[[Linc_RNA]] )){
     ProCod_pos <- "NA"
     if(is.null(proNeg_hash[[Linc_RNA]])){
          ProCod_neg <- "NA"
          GenoCordinates <- "0000"
          linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
        }
      else{
          ProCod_neg <- proNeg_hash[[Linc_RNA]]
          if(Start_linc < Start_negPro){
            GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_linc),"-",as.character(End_negPro))
            linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
            }
          else if ( Start_negPro < Start_linc ){
            GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_negPro),"-",as.character(End_linc))
            linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
            }
        }
    }
  else {
    ProCod_pos <- proPos_hash[[Linc_RNA]]
    if(is.null(proNeg_hash[[Linc_RNA]])){
      ProCod_neg <- "NA"
        if(Start_linc < Start_posPro){
        GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_linc),"-",as.character(End_posPro))
        linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
        }
        else if ( Start_posPro < Start_linc ){
          GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_posPro),"-",as.character(End_linc))
          linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
          }
      }
    else{
      ProCod_neg <- proNeg_hash[[Linc_RNA]]
      if((Start_posPro < Start_linc) && (Start_negPro < Start_linc)){
        if(Start_posPro < Start_negPro){
          GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_posPro),"-",as.character(End_linc))
          linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
        }
        if(Start_negPro < Start_posPro){
          GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_negPro),"-",as.character(End_linc))
          linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
          }
      }
      if((Start_posPro > Start_linc) && (Start_negPro > Start_linc)){
        if(Start_posPro < Start_negPro){
          GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_linc),"-",as.character(End_negPro))
          linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
          }
        if(Start_negPro < Start_posPro){
            GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_linc),"-",as.character(End_posPro))
            linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
          }
      }
      if((Start_posPro > Start_linc) && (Start_negPro < Start_linc)){
          GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_negPro),"-",as.character(End_posPro))
          linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
          }
      if((Start_posPro < Start_linc) && (Start_negPro > Start_linc)){
          GenoCordinates <- paste0("chr",Chrom,":",as.character(Start_posPro),"-",as.character(End_negPro))
          linc_Pro_pairs <- rbind(linc_Pro_pairs, data.frame(Linc_RNA, ProCod_pos, ProCod_neg, GenoCordinates))
        }
      }
    }
}
write.table(linc_Pro_pairs, "/home/layal/OSTEO_2017/lnc_RNA/lincs_proPos_proNeg.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# ============================================================ Heatmap ==========================================================================================
allLincs <- rownames(DifGE_linc)
DSC <- as.data.frame(assay(ddseq)[allLincs,])
# mean of replivcates
DSCM <-data.frame(row.names = rownames(DSC))
DSCM$Day0_Mean <- rowMeans(cbind(DSC$DayZero_R1, DSC$DayZero_R2 , DSC$DayZero_R3), na.rm=TRUE)
DSCM$Day3_Mean <- rowMeans(cbind(DSC$DayThree_R1 , DSC$DayThree_R2 , DSC$DayThree_R3) , na.rm = TRUE)
DSCM$Day6_Mean <- rowMeans(cbind(DSC$DaySix_R1 , DSC$DaySix_R2, DSC$DaySix_R3) , na.rm = TRUE)
DSCM$Day12_Mean <- rowMeans(cbind(DSC$DayTwelve_R1, DSC$DayTwelve_R2 , DSC$DayTwelve_R3) , na.rm = TRUE)
# scaling new value: (x – μ)/σ where x is the value , μ is the mean of the rows, σ is the stadard diviation, however the mean of the new row values is 0 
DSCM_scale <- apply(DSCM, 1, scale)

# we need to transpose the scales matrix so we get the genes for clustering
myclusters <-t(DSCM_scale)
colnames(myclusters) <- colnames(DSCM)

pdf(paste(path,"Heatmap-MeanLinc.pdf",sep = "/"),height=60,width=20)
heatmap.2(myclusters,dendrogram="row",col=rev(rainbow(300, start=0, end=4/6)),margins = c(25,10), 
          Colv=F,Rowv=T,trace="none",key=T,cexCol = 4,cexRow = 0.9, scale = "none", keysize=0.7)
abline(v = c(0.32,0.525,0.73), untf = FALSE, col = "black",lwd=3)
abline(v = 0.005, untf = FALSE, col = "red")
dev.off()

#====================================== Getting the clusters ======================================================

#----------------------------------------cluster A--------------------------------------
cluster_A <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/A.txt", header = F)
colnames(cluster_A) <- "LincRNA"

A_pos <- posCor[which(posCor$LincRNA %in% cluster_A$LincRNA), ]
A_posPro <- as.character(A_pos$Protein_cod)

A_neg <- negCor[which(negCor$LincRNA %in% cluster_A$LincRNA) , ]
A_negPro <-as.character( A_neg$Protein_cod)

A_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_A$LincRNA) ,]
A_notAnnPro <- c(rep("NA", length(rownames(A_notAnn))))

write(A_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/A_posPro.txt")
write(A_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/A_negPro.txt")
write(A_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/A_notAnnPro.txt")

#-----------------------------------cluster B----------------------------------------------------
cluster_B <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/B.txt", header = F)
colnames(cluster_B) <- "LincRNA"

B_pos <- posCor[which(posCor$LincRNA %in% cluster_B$LincRNA), ]
B_posPro <-as.character( B_pos$Protein_cod)

B_neg <- negCor[which(negCor$LincRNA %in% cluster_B$LincRNA) , ]
B_negPro <- as.character( B_neg$Protein_cod)
B_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_B$LincRNA) ,]
B_notAnnPro <- c(rep("NA", length(rownames(B_notAnn))))

write(B_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/B_posPro.txt")
write(B_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/B_negPro.txt")
write(B_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/B_notAnnPro.txt")

#-----------------------------------cluster C----------------------------------------------------
cluster_C <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/C.txt", header = F)
colnames(cluster_C) <- "LincRNA"

C_pos <- posCor[which(posCor$LincRNA %in% cluster_C$LincRNA), ]
C_posPro <- as.character(C_pos$Protein_cod)

C_neg <- negCor[which(negCor$LincRNA %in% cluster_C$LincRNA) , ]
C_negPro <- as.character( C_neg$Protein_cod)

C_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_C$LincRNA) ,]
C_notAnnPro <- c(rep("NA", length(rownames(C_notAnn))))

write(C_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/C_posPro.txt")
write(C_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/C_negPro.txt")
write(C_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/C_notAnnPro.txt")


#-----------------------------------cluster D----------------------------------------------------
cluster_D <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/D.txt", header = F)
colnames(cluster_D) <- "LincRNA"

D_pos <- posCor[which(posCor$LincRNA %in% cluster_D$LincRNA), ]
D_posPro <- as.character( D_pos$Protein_cod)

D_neg <- negCor[which(negCor$LincRNA %in% cluster_D$LincRNA) , ]
D_negPro <- as.character(D_neg$Protein_cod)
D_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_D$LincRNA) ,]
D_notAnnPro <- c(rep("NA", length(rownames(D_notAnn))))

write(D_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/D_posPro.txt")
write(D_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/D_negPro.txt")
write(D_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/D_notAnnPro.txt")


#-----------------------------------cluster E----------------------------------------------------
cluster_E <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/E.txt", header = F)
colnames(cluster_E) <- "LincRNA"

E_pos <- posCor[which(posCor$LincRNA %in% cluster_E$LincRNA), ]
E_posPro <- as.character(E_pos$Protein_cod)

E_neg <- negCor[which(negCor$LincRNA %in% cluster_E$LincRNA) , ]
E_negPro <- as.character(E_neg$Protein_cod)
E_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_E$LincRNA) ,]
E_notAnnPro <- c(rep("NA", length(rownames(E_notAnn))))

write(E_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/E_posPro.txt")
write(E_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/E_negPro.txt")
write(E_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/E_notAnnPro.txt")

#-----------------------------------cluster F----------------------------------------------------
cluster_F <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/F.txt", header = F)
colnames(cluster_F) <- "LincRNA"

F_pos <- posCor[which(posCor$LincRNA %in% cluster_F$LincRNA), ]
F_posPro <- as.character(F_pos$Protein_cod)

F_neg <- negCor[which(negCor$LincRNA %in% cluster_F$LincRNA) , ]
F_negPro <- as.character(F_neg$Protein_cod)
F_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_F$LincRNA) ,]
F_notAnnPro <- c(rep("NA", length(rownames(F_notAnn))))

write(F_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/F_posPro.txt")
write(F_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/F_negPro.txt")
write(F_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/F_notAnnPro.txt")

#-----------------------------------cluster G----------------------------------------------------
cluster_G <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/G.txt", header = F)
colnames(cluster_G) <- "LincRNA"

G_pos <- posCor[which(posCor$LincRNA %in% cluster_G$LincRNA), ]
G_posPro <- as.character(G_pos$Protein_cod)

G_neg <- negCor[which(negCor$LincRNA %in% cluster_G$LincRNA) , ]
G_negPro <- as.character(G_neg$Protein_cod)
G_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_G$LincRNA) ,]
G_notAnnPro <- c(rep("NA", length(rownames(G_notAnn))))

write(G_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/G_posPro.txt")
write(G_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/G_negPro.txt")
write(G_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/G_notAnnPro.txt")

#-----------------------------------cluster H----------------------------------------------------
cluster_H <- read.csv("/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/H.txt", header = F)
colnames(cluster_H) <- "LincRNA"

H_pos <- posCor[which(posCor$LincRNA %in% cluster_H$LincRNA), ]
H_posPro <- as.character(H_pos$Protein_cod)

H_neg <- negCor[which(negCor$LincRNA %in% cluster_H$LincRNA) , ]
H_negPro <- as.character( H_neg$Protein_cod)
H_notAnn <- notLinc[which(rownames(notLinc) %in% cluster_H$LincRNA) ,]
H_notAnnPro <- c(rep("NA", length(rownames(H_notAnn))))

write(H_posPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/H_posPro.txt")
write(H_negPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/H_negPro.txt")
write(H_notAnnPro, "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/lnc_RNA/GO_linc/GO_pro/H_notAnnPro.txt")



