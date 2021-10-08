library("ggplot2")
library("reshape2")
library("grid")
library("gridExtra")
install.packages("cowplot")
library("cowplot")
library("DESeq2")

path="/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq"
load(paste(path,"DESeq.RData", sep = "/"))


x <- data.frame( Day12_0 = res12_0$log2FoldChange , Day12_3 = res12_3$log2FoldChange , Day12_6 = res12_6$log2FoldChange ,
                 Day6_0 = res6_0$log2FoldChange ,Day3_0 = res3_0$log2FoldChange ,Day6_3 = res6_3$log2FoldChange )
head(x)
dataDensity<- melt(x)
head(dataDensity)

x12 <- data.frame( Day12_0 = res12_0$log2FoldChange )
dataDensity12<- melt(x12)
p12 <- ggplot(dataDensity12 , aes(x =value)) + geom_density(alpha=0.75, fill="#0000FF") +ggtitle("Day12 vs Day0")+xlim(-4,4)+ylim(0.0,0.7)+ylab("")+xlab("")+theme(plot.title = element_text(size=12, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=12))

x6 <- data.frame( Day6_0 = res6_0$log2FoldChange )
dataDensity6<- melt(x6)
p6 <- ggplot(dataDensity6 , aes(x =value)) + geom_density(alpha=0.75, fill="#0000FF")+ggtitle("Day6 vs Day0") +xlim(-4,4)+ylim(0.0,0.7)+ylab("")+xlab("")+theme(plot.title = element_text(size=12, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=12))

x3 <- data.frame( Day3_0 = res3_0$log2FoldChange )
dataDensity3<- melt(x3)
p3 <- ggplot(dataDensity3 , aes(x =value)) + geom_density(alpha=0.75, fill="#0000FF") +ggtitle("Day3 vs Day0")+xlim(-4,4)+ylim(0.0,0.8)+ylab("")+xlab("")+theme(plot.title = element_text(size=12, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=12))

x12_3 <- data.frame( Day12_3 = res12_3$log2FoldChange )
dataDensity12_3<- melt(x12_3)
p12_3 <- ggplot(dataDensity12_3 , aes(x =value)) + geom_density(alpha=0.75, fill="#0000FF") +ggtitle("Day12 vs Day3")+xlim(-4,4)+ylim(0.0,0.9)+ylab("")+xlab("")+theme(plot.title = element_text(size=12, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=12))

x12_6 <- data.frame( Day12_6 = res12_6$log2FoldChange )
dataDensity12_6<- melt(x12_6)
p12_6 <- ggplot(dataDensity12_6 , aes(x =value)) + geom_density(alpha=0.75, fill="#0000FF") +ggtitle("Day12 vs Day6")+xlim(-4,4)+ylim(0.0,0.95)+ylab("")+xlab("")+theme(plot.title = element_text(size=12, face = "bold"),axis.text=element_text(size=12), axis.title=element_text(size=12))

x6_3 <- data.frame( Day6_3 = res6_3$log2FoldChange )
dataDensity6_3<- melt(x6_3)
p6_3 <- ggplot(dataDensity6_3 , aes(x =value)) + geom_density(alpha=0.75, fill="#0000FF")+ggtitle("Day6 vs Day3")+xlim(-3,3)+ylim(0.0,1.6)+ylab("")+xlab("")+theme(plot.title = element_text(size=12, face = "bold"), axis.text=element_text(size=12), axis.title=element_text(size=12))


mgp<- grid.arrange(p3,p6,p12,p6_3,p12_3,p12_6, ncol=3, nrow =2,
                   top=textGrob("Fold Change distribution across all comparisons", gp=gpar(fontface="bold",fontsize=16)), 
                   left=textGrob("Density", rot=90 ,gp=gpar(fontface="bold",fontsize=16)),
                   bottom= textGrob("log2(FC)" ,gp=gpar(fontface="bold",fontsize=16)))

ggsave(paste(path,"multiFCdensity.pdf", sep = "/"), mgp, scale = 1, width=18, height=16, units = "cm")
