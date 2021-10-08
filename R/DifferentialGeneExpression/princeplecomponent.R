library(DESeq2)
library(ggplot2)

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

#===============================  regularized logarithm =============================================
rld <- rlog(ddseq, blind = FALSE)
dataPCA_rlg <- plotPCA(rld, intgroup= "condition", ntop=length(rownames(ddseq)), returnData=TRUE)
percentVar <- round(100 * attr(dataPCA_rlg, "percentVar"))

pdf( "PCA_ggplot_rlg.pdf", width = 3.35, height = 3.5, onefile=FALSE)
ggplot(dataPCA_rlg, aes(PC1, PC2, shape=condition)) +
  geom_point(size=4)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_shape_discrete(solid=T, legend=T) +
  ggtitle("Principal Component") + 
  theme(plot.title = element_text(face = "bold",size = 11, hjust = 0.5, family = 'sans'),
        text = element_text(size = 10),
        panel.background = element_blank(),
        axis.title = element_text(face = "bold",size = 10, family = 'sans'),
        axis.text.x = element_text(size = 10, colour = 'black', family = 'sans'),
        axis.text.y = element_text(size = 10, colour = 'black', family = 'sans'),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        legend.key = element_rect(colour = 'white'),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(face= "italic" ,size = 10, colour = 'black', family = 'sans'),
        legend.title = element_blank(),
        plot.margin = margin(0.3,0.2,0.05,0.3,"cm")
  )
dev.off()

#============================ variance Stabilizing Transformation =============================================
vsdd <- varianceStabilizingTransformation(ddseq , blind = FALSE)
dataPCA <- plotPCA(vsdd, intgroup="condition", ntop=length(rownames(ddseq)), returnData=TRUE)
percentVar1 <- round(100 * attr(dataPCA, "percentVar"))

pdf( "PCA_ggplot_vst.pdf", width = 3.35, height = 3.5, onefile=FALSE)
ggplot(dataPCA, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) + 
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
  ggtitle("Principal Component") + 
  theme(plot.title = element_text(face = "bold",size = 11, hjust = 0.5, family = 'sans'),
        text = element_text(size = 10),
        panel.background = element_blank(),
        axis.title = element_text(face = "bold",size = 10, family = 'sans'),
        axis.text.x = element_text(size = 10, colour = 'black', family = 'sans'),
        axis.text.y = element_text(size = 10, colour = 'black', family = 'sans'),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        legend.key = element_rect(colour = 'white'),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(face= "italic" ,size = 10, colour = 'black', family = 'sans'),
        legend.title = element_blank(),
        plot.margin = margin(0.3,0.3,0.05,0.3,"cm")
  )


dev.off()
