library (ggplot2)
library(DESeq2)

path = "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq"
load(paste(path, "DESeq.RData", sep = "/"))

pdf(file = paste(path,"volcanoAllplots.pdf",sep = "/"), height=5.5,width=6.7)
my_FC <- log2(5)
my_padj <- 10^(-50)

layout(matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3,ncol = 3, byrow=TRUE), heights=c(4,4, 1))
par(pch = 20 ,cex=0.5,cex.axis=1.5, cex.main = 1.8, oma = c(2, 5, 0.5, 0.5), mar=c(2,1.6,2.5,1.5)) 

#---------------p1-----------
with(resThree_Zero, plot(log2FoldChange, -log10(padj), main = "Day3 vs Day0", col = "blue", ylim=c(0,300) , xlim = c(-8,8) ,
                         xlab="",ylab ="" ))
with(subset(resThree_Zero, padj < my_padj & abs(log2FoldChange) > my_FC), points(log2FoldChange,-log10(padj), col = "red"))
#---------------p2-----------
with(resSix_Zero, plot(log2FoldChange, -log10(padj), main = "Day6 vs Day0", col = "blue", ylim=c(0,300) , xlim = c(-8,8),
                       xlab="",ylab =""))
with(subset(resSix_Zero, padj < my_padj & abs(log2FoldChange) > my_FC), points(log2FoldChange,-log10(padj), col = "red"))
#---------------p3-----------
with(resTwelve_Zero, plot(log2FoldChange, -log10(padj), main = "Day12 vs Day0", col = "blue", ylim=c(0,300) , xlim = c(-8,8),
                          xlab="",ylab =""))
with(subset(resTwelve_Zero, padj < my_padj & abs(log2FoldChange) > my_FC), points(log2FoldChange,-log10(padj), col = "red"))
#---------------p4-----------
with(resTwelve_Three, plot(log2FoldChange, -log10(padj), main = "Day12 vs Day3", col = "blue", ylim=c(0,300) , xlim = c(-8,8),
                           xlab="",ylab =""))
with(subset(resTwelve_Three, padj < my_padj& abs(log2FoldChange) > my_FC), points(log2FoldChange,-log10(padj), col = "red"))
#---------------p5-----------
with(resTwelve_Six, plot(log2FoldChange, -log10(padj), main = "Day12 vs Day6", col = "blue",ylim=c(0,300), xlim = c(-8,8),
                         xlab="",ylab =""))
with(subset(resTwelve_Six, padj < my_padj & abs(log2FoldChange) > my_FC), points(log2FoldChange,-log10(padj), col = "red"))
#---------------p6-----------
with(resSix_Three, plot(log2FoldChange, -log10(padj), main = "Day6 vs Day3", col = "blue",ylim=c(0,300) , xlim = c(-8,8),
                        xlab="",ylab =""))
with(subset(resSix_Three, padj < my_padj & abs(log2FoldChange) > my_FC), points(log2FoldChange,-log10(padj), col = "red"))

par(mai=c(0,0,0,0))
plot.new()
legend(x="left", ncol=1,legend = c("All DE Genes", "(adj_P < 10^-50) & (FC > 5)" ),
       col =c("blue","red"), cex=1.6 ,pt.cex = 2, pch=16)
# common axis
mtext(text="log2 Fold Change",side=1,line=-3,outer=TRUE, cex = 1.3 )
mtext(text="-log10 adjusted P_value ",side=2,line=2,outer=TRUE, cex = 1.3)

dev.off()

#================================================================================

pdf(file = "volcano3_0.pdf", height=6.5,width=6.2)
FC <- log2(4)
pj <- 0.001
layout(matrix(c(1,2,3,4,5,5), nrow = 3,ncol = 2, byrow=TRUE), heights=c(4,4, 1.6))
par(pch = 20 ,cex=0.5,cex.axis=1.5, cex.main = 1.8, oma = c(2, 5, 0.5, 0.5), mar=c(2,1.6,2.5,1.5)) 

with(res3_0, plot(log2FoldChange, -log10(padj), col = "black", ylim=c(0,300) , xlim = c(-8,8) ,
                         xlab="",ylab ="" ))
with(subset(res3_0, padj < pj), points(log2FoldChange,-log10(padj), col = "blue"))
title("(A)", adj = 0, line = 1,cex=1.2)

with(resThree_Zero, plot(log2FoldChange, -log10(padj), col = "black", ylim=c(0,300) , xlim = c(-8,8) ,
                         xlab="",ylab ="" ))
with(subset(resThree_Zero, padj < pj), points(log2FoldChange,-log10(padj), col = "blue"))
with(subset(resThree_Zero, padj < pj & abs(log2FoldChange) > FC), points(log2FoldChange,-log10(padj), col = "red"))
title("(B)", adj = 0, line = 1,cex=1.2)

with(resThree_Zero, plot(log2FoldChange, -log10(padj), col = "black", ylim=c(0,300) , xlim = c(-8,8) ,
                         xlab="",ylab ="" ))
with(subset(resThree_Zero, padj < my_padj), points(log2FoldChange,-log10(padj), col = "blue"))
title("(C)", adj = 0, line = 1,cex=1.2)

with(resThree_Zero, plot(log2FoldChange, -log10(padj), col = "black", ylim=c(0,300) , xlim = c(-8,8) ,
                         xlab="",ylab ="" ))
with(subset(resThree_Zero, padj < my_padj), points(log2FoldChange,-log10(padj), col = "blue"))
with(subset(resThree_Zero, padj < my_padj & abs(log2FoldChange) > my_FC), points(log2FoldChange,-log10(padj), col = "red"))
title("(D)", adj = 0, line = 1,cex=1.2)

plot.new()
legend(x="left", ncol=1,legend = c("All DE Genes", "p-value cutoff","p-value & LFC cutoffs" ),
       col =c("black","blue","red"), cex=1.6 ,pt.cex = 2, pch=16)
# common axis
mtext(text="log2 Fold Change",side=1,line=-8,outer=TRUE, cex = 1.3 )
mtext(text="-log10 adjusted P_value ",side=2,line=2,outer=TRUE, cex = 1.3)
dev.off()


