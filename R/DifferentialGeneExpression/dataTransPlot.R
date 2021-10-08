library(DESeq2)
library(ggplot2)
library("vsn")
library(grid)
library(gridExtra)
#Note that the vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions. While a flat curve of the square root of variance over the mean may seem like the goal of such transformations, this may be unreasonable in the case of datasets with many true differences due to the experimental conditions.

pathOut = "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq"

FC_countreads <- readRDS(paste(pathOut,"FC_counts.rds",sep = "/"))
countrawdata <- FC_countreads$counts
countdata <- countrawdata
colnames(countdata) <- gsub(pattern = "X.home.layal.Solexa.layal.OSTEO.Align_Tophat_GRCm38.|.accepted_hits.bam",replacement = "",x = colnames(countdata))
head(countdata)
New_names <- c("Day0_R1","Day0_R2","Day0_R3" ,"Day3_R1","Day3_R2","Day3_R3","Day6_R1","Day6_R2","Day6_R3", "Day12_R1", "Day12_R2", "Day12_R3")
colnames(countdata) <- New_names
condition <- factor(c(rep("Day0", 3),rep("Day3",3),rep("Day6",3),rep("Day12",3)), levels = c("Day0", "Day3","Day6", "Day12"))
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
ddsf <- dds[ rowSums(counts(dds)) > 1, ]
ddseq <- DESeq(ddsf)

ntd <- normTransform(ddseq)
rld <- rlog(ddseq, blind=FALSE)
vsd <- varianceStabilizingTransformation(ddseq, blind=FALSE)
notAllZero <- (rowSums(counts(ddseq))>0)

gnt <- meanSdPlot(assay(ntd)[notAllZero,] , show.legend = FALSE, xlab  = "", ylab  = "" )
grl <- meanSdPlot(assay(rld)[notAllZero,] , show.legend = FALSE, xlab  = "", ylab  = "" )
gvs <- meanSdPlot(assay(vsd)[notAllZero,] , show.legend = FALSE, xlab  = "", ylab  = "" )

themePlot <- function(base_size=12, base_family="sans") {
    theme(plot.title = element_text(face = "bold",size = 12, hjust = 0, family = 'sans'),
        text = element_text(size = 11),
        axis.text.x = element_text(size = 12, colour = 'black', family = 'sans'),
        axis.text.y = element_text(size = 12, colour = 'black', family = 'sans'),
        axis.line = element_line(colour="black", size = 0.5),
        legend.text = element_text(face= "italic" ,size = 12, colour = 'black', family = 'sans'))
}


g1 <- gnt$gg + ggtitle("shifted logarithm transformation")+themePlot()
g2 <- grl$gg + ggtitle("regularized log transformation")+themePlot()
g3 <- gvs$gg + ggtitle("variance stabilizing transformation")+themePlot()

grid.arrange(g1,g2, g3, ncol=3)
grid.text("Mean", x=0.5, y=0.03, gp=gpar(fontsize=16))
grid.text("SD", x=0.01, y=0.6, rot=90, gp=gpar(fontsize=16))
