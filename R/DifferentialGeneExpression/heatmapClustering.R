# biocLite("DESeq2")
#citation("DESeq2")
# biocLite("pheatmap")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(gplots)

path = "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq"
load(paste(path, "DESeq.RData", sep = "/"))

genesCluster <- genesProt 

#================================== heatmap of transformation counts  =====================
vsd <- varianceStabilizingTransformation(ddseq, blind=FALSE)
DSV <- assay(vsd)[genesCluster,]
myclustersVT <- t(apply(as.matrix(DSV), 1, scale))
colnames(myclustersVT) <- colnames(DSV)

pdf(paste(path,"HeatmapVST.pdf",sep = "/"),height=60,width=20)
heatmap.2(myclustersVT,dendrogram="row",col=rev(rainbow(300, start=0, end=4/6)),margins = c(25,10), 
          Colv=F,Rowv=T,trace="none",key=T,cexCol = 4,cexRow = 0.5, scale = "none", keysize=0.7)
abline(v = c(0.32,0.525,0.73), untf = FALSE, col = "black",lwd=3)
abline(v = 0.01, untf = FALSE, col = "red")
dev.off()


#--------------------Mean of replicates
DSCV <- as.data.frame(assay(vsd)[genesCluster,])
# mean of replivcates
DSCMv <-data.frame(row.names = rownames(DSCV))
DSCMv$Day0_Mean <- rowMeans(cbind(DSCV$DayZero_R1, DSCV$DayZero_R2 , DSCV$DayZero_R3), na.rm=TRUE)
DSCMv$Day3_Mean <- rowMeans(cbind(DSCV$DayThree_R1 , DSCV$DayThree_R2 , DSCV$DayThree_R3) , na.rm = TRUE)
DSCMv$Day6_Mean <- rowMeans(cbind(DSCV$DaySix_R1 , DSCV$DaySix_R2, DSCV$DaySix_R3) , na.rm = TRUE)
DSCMv$Day12_Mean <- rowMeans(cbind(DSCV$DayTwelve_R1, DSCV$DayTwelve_R2 , DSCV$DayTwelve_R3) , na.rm = TRUE)
# scaling new value: (x – μ)/σ where x is the value , μ is the mean of the rows, σ is the stadard diviation, however the mean of the new row values is 0 
DSCMv_scale <- apply(DSCMv, 1, scale)

# we need to transpose the scales matrix so we get the genes for clustering
myclustersVM <-t(DSCMv_scale)
colnames(myclustersVM) <- colnames(DSCMv)

pdf(paste(path,"Heatmap-MeanVST.pdf",sep = "/"),height=60,width=20)
heatmap.2(myclustersVM,dendrogram="row",col=rev(rainbow(300, start=0, end=4/6)),margins = c(25,10), 
          Colv=F,Rowv=T,trace="none",key=T,cexCol = 4,cexRow = 0.9, scale = "none", keysize=0.7)
abline(v = c(0.32,0.525,0.73), untf = FALSE, col = "black",lwd=3)
abline(v = 0.01, untf = FALSE, col = "red")
dev.off()

# =============================== Heatmap of ddseq with replicates =======================================================

# ---------------Calculate the Variance 

DSC <- assay(ddseq)[genesCluster,]
DFgenes_var = apply(DSC, 1, var)

size = length(DFgenes_var)
var_select = order(DFgenes_var, decreasing = T)[1:size]
myhigh.var.genes = DSC[var_select, ]
myclusters <- t(apply(as.matrix(myhigh.var.genes), 1, scale))
colnames(myclusters) <- colnames(DSC)

pdf(paste(path,"HeatmapAllFC5.pdf",sep = "/"),height=60,width=20)
heatmap.2(myclusters,dendrogram="row",col=rev(rainbow(300, start=0, end=5/6)),margins = c(25,10), 
          Colv=F,Rowv=T,trace="none",key=T,cexCol = 4,cexRow = 1, scale = "none", keysize=0.7)
dev.off()

# ------------- without calculating variance --------------------
DSC <- assay(ddseq)[genesCluster,]
myclusters <- t(apply(as.matrix(DSC), 1, scale))
colnames(myclusters) <- colnames(DSC)

pdf(paste(path,"HeatmapAll2.pdf",sep = "/"),height=60,width=20)
heatmap.2(myclusters,dendrogram="row",col=rev(rainbow(300, start=0, end=5/6)),margins = c(25,10), 
          Colv=F,Rowv=T,trace="none",key=T,cexCol = 4,cexRow = 0.5, scale = "none", keysize=0.7)
abline(v = c(0.32,0.525,0.73), untf = FALSE, col = "black",lwd=3)
abline(v = 0.01, untf = FALSE, col = "red")
dev.off()

NCOF <- as.data.frame(NormSF_Counts[which(rownames(NormSF_Counts) %in% OnOffGenes), ]) 
COFG <-data.frame(row.names = rownames(NCOF))
#--------------------Mean of replicates
DSC <- as.data.frame(assay(ddseq)[genesCluster,])
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

pdf(paste(path,"Heatmap-Mean.pdf",sep = "/"),height=60,width=20)
heatmap.2(myclusters,dendrogram="row",col=rev(rainbow(300, start=0, end=4/6)),margins = c(25,10), 
          Colv=F,Rowv=T,trace="none",key=T,cexCol = 4,cexRow = 0.9, scale = "none", keysize=0.7)
abline(v = c(0.32,0.525,0.73), untf = FALSE, col = "black",lwd=3)
abline(v = 0.01, untf = FALSE, col = "red")
dev.off()






# ============================================== Hierarchical clustering analysis ============================================================

## Get the sample distances from the transpose of the whole thing
rld <- rlog(ddseq, blind=FALSE)
sampledist <- dist(t(assay(rld)))
pdf(paste(path,"Hierarchical_clustering.pdf",sep = "/"))
plot(hclust(sampledist))
dev.off()

library("RColorBrewer")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampledistmat <- as.matrix(sampledist)
hc <- hclust(sampledist)
pdf(paste(path,"Heatmap_distancesBlue.pdf",sep = "/"))
heatmap.2(sampledistmat,col=colors, key = TRUE, trace = "none",margins = c(7, 7), keysize=1.5 , dendrogram="both")
dev.off()

pdf(paste(path,"Heatmap_distancesBlue_order.pdf",sep = "/"))
heatmap.2(sampledistmat,col=colors, key = TRUE, trace = "none",margins = c(8, 8), keysize=1.5, Rowv = "as-is" ,Colv ="RowV" , dendrogram="none")
dev.off()

pdf(paste(path,"Heatmap_distancesBlue_order_Dendo.pdf",sep = "/"))
heatmap.2(sampledistmat,col=colors, key = TRUE, trace = "none",margins = c(7, 7), keysize=1.5, Rowv = "as-is" , dendrogram="both")
dev.off()
