# biocLite() is the bioconductor installer function. Run it without any
# arguments to install the core packages or update any installed packages.
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("DESeq2")
library(DESeq2)
citation("DESeq2")
library(ggplot2)

pathMA="/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/DESeq"
load(paste(path, "DESeq.RData", sep = "/"))


resTwelve_Zero <- results(ddseq, contrast=c("condition", "DayTwelve", "DayZero"),pAdjustMethod = "bonferroni", alpha = 0.01,
                          lfcThreshold=2, altHypothesis="greaterAbs")
resSix_Zero <- results(ddseq, contrast=c("condition", "DaySix" ,"DayZero") , pAdjustMethod = "bonferroni", alpha = 0.01,
                       lfcThreshold=2, altHypothesis="greaterAbs")
resThree_Zero <- results(ddseq, contrast=c("condition", "DayThree", "DayZero"), pAdjustMethod = "bonferroni", alpha = 0.01,
                         lfcThreshold=2, altHypothesis="greaterAbs")
resTwelve_Three <- results(ddseq, contrast=c("condition", "DayTwelve", "DayThree"),pAdjustMethod = "bonferroni", alpha = 0.01,
                           lfcThreshold=2, altHypothesis="greaterAbs")
resTwelve_Six <- results(ddseq, contrast=c("condition", "DayTwelve", "DaySix"),pAdjustMethod = "bonferroni", alpha = 0.01,
                         lfcThreshold=2, altHypothesis="greaterAbs")
resSix_Three <- results(ddseq, contrast=c("condition", "DaySix", "DayThree"),pAdjustMethod = "bonferroni", alpha = 0.01,
                        lfcThreshold=2, altHypothesis="greaterAbs")


pdf(file = paste(pathMA,"MA_plots0.01_LFC2.pdf",sep = "/"), height=5.5,width=6.7)
a =0.01
layout(matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3,ncol = 3, byrow=TRUE), heights=c(4.1,4.1, 1))
par(pch = 20 ,cex=0.6, cex.axis=1.115, cex.main = 1.7, oma = c(1, 4.5, 4, 0.5), mar=c(2,1.5,3,1.5))

plotMA(resThree_Zero, alpha= a, main="D3_D0", ylim=c(-6,6), xlab= "", ylab="")
abline(h=c(-2,2),col="dodgerblue",lwd=2)

plotMA(resSix_Zero, alpha= a, main="D6_D0", ylim=c(-6,6), xlab= "", ylab="")
abline(h=c(-2,2),col="dodgerblue",lwd=2)

plotMA(resTwelve_Zero, alpha= a,  main=" D12_D0", ylim=c(-6,6), xlab= "", ylab="")
abline(h=c(-2,2),col="dodgerblue",lwd=2)

plotMA(resTwelve_Three, alpha= a , main=" D12_D3", ylim=c(-6,6), xlab= "", ylab="")
abline(h=c(-2,2),col="dodgerblue",lwd=2)

plotMA(resTwelve_Six, alpha= a, main="D12_D6", ylim=c(-6,6), xlab= "", ylab="")
abline(h=c(-2,2),col="dodgerblue",lwd=2)

plotMA(resSix_Three, alpha= a, main="D6_D3", ylim=c(-6,6), xlab= "", ylab="")
abline(h=c(-2,2),col="dodgerblue",lwd=2)
par(mai=c(0,0,0,0))
plot.new()
legend(x="left", ncol=1,legend = "alpha =< 0.01" ,col ="red", cex=1.4 ,pt.cex = 2.5, pch=16)

# common axis
mtext(text="Mean of Normalized Counts",side=1,line=-3,outer=TRUE, cex = 1.2 )
mtext(text=" Fold change (log2) ",side=2,line=2,outer=TRUE, cex = 1.2)

dev.off()

