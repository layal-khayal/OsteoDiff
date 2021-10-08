library(DESeq2)
library(ggplot2)
install.packages("pheatmap")
library(pheatmap)
library(gplots)

pathOut = "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq"
FC_countreads <- readRDS("FC_counts.rds")
countrawdata <- FC_countreads$counts
countdata <- countrawdata
colnames(countdata) <- gsub(pattern = "X.home.layal.Solexa.layal.OSTEO.Align_Tophat_GRCm38.|.accepted_hits.bam",replacement = "",x = colnames(countdata))
New_names <- c("Day0_R1","Day0_R2","Day0_R3" ,"Day3_R1","Day3_R2","Day3_R3","Day6_R1","Day6_R2","Day6_R3", "Day12_R1", "Day12_R2", "Day12_R3")
colnames(countdata) <- New_names

condition <- factor(c(rep("Day0", 3),rep("Day3",3),rep("Day6",3),rep("Day12",3)), levels = c("Day0", "Day3","Day6", "Day12"))
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
ddsf <- dds[ rowSums(counts(dds)) > 1, ]
ddseq <- DESeq(ddsf)

## Get the sample distances from the transpose of the whole thing
rld <- rlog(ddseq, blind=FALSE)
sampledist <- dist(t(assay(rld)))


pdf(paste(pathOut,"Hierarchical_clustering.pdf",sep ='/'), width = 4, height = 3.5, onefile=FALSE )
x <- hclust(sampledist)
par(cex=0.7, cex.axis=1, cex.lab=1, cex.main = 1.2, oma = c(1,1,1,1), mar=c(0,4,2,2))
plot(x, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Dendrogram of Osteoblast dataset",ylab = "Height",
     sub = "", xlab ="")
dev.off()

#library("RColorBrewer")
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
n=20
gray.colors(n, start = 0, end = 1)
sampledistmat <- as.matrix(sampledist)
colnames(sampledistmat) <- NULL
pdf(paste(path,"Samples_distancesBlue.pdf", sep = '/'), width = 4.5, height = 3.5, onefile=FALSE)
pheatmap(sampledistmat,
         clustering_distance_rows=sampledist,
         clustering_distance_cols=sampledist,
         col=colors,
         fontsize = 8)
dev.off()


#heatmap.2(sampledistmat,cexRow = 0.6,cexCol=0.6, col=colors, key = FALSE, trace = "none", symm =TRUE, margins = c(3.5,3.5), dendrogram='column')