library(DESeq2)
library(ggplot2)
path="/home/layal/OSTEO_2017/DESeq"
load(paste(path, "DESeq.RData", sep = "/"))

pathClus <- "/home/layal/OSTEO_2017/GeneOntology"

clustA <- scan(paste(pathClus,"ClusterA.txt",sep="/"), what = character(),sep = "\n")
clustB <- scan(paste(pathClus,"ClusterB.txt",sep="/"), what = character(),sep = "\n")
clustC <- scan(paste(pathClus,"ClusterC.txt",sep="/"), what = character(),sep = "\n") #CD1
clustD <- scan(paste(pathClus,"ClusterD.txt",sep="/"), what = character(),sep = "\n") #D2
clustE <- scan(paste(pathClus,"ClusterE.txt",sep="/"), what = character(),sep = "\n")
clustF <- scan(paste(pathClus,"ClusterF.txt",sep="/"), what = character(),sep = "\n")
clustG <- scan(paste(pathClus,"ClusterG.txt",sep="/"), what = character(),sep = "\n")
clustH <- scan(paste(pathClus,"ClusterH.txt",sep="/"), what = character(),sep = "\n")
clustI <- scan(paste(pathClus,"ClusterI.txt",sep="/"), what = character(),sep = "\n")

genesCluster <- genesProt 
DSC <- as.data.frame(assay(ddseq)[genesCluster,])
# means of replicates
DSCM <-data.frame(row.names = rownames(DSC))
DSCM$Day0_Mean <- rowMeans(cbind(DSC$DayZero_R1, DSC$DayZero_R2 , DSC$DayZero_R3), na.rm=TRUE)
DSCM$Day3_Mean <- rowMeans(cbind(DSC$DayThree_R1 , DSC$DayThree_R2 , DSC$DayThree_R3) , na.rm = TRUE)
DSCM$Day6_Mean <- rowMeans(cbind(DSC$DaySix_R1 , DSC$DaySix_R2, DSC$DaySix_R3) , na.rm = TRUE)
DSCM$Day12_Mean <- rowMeans(cbind(DSC$DayTwelve_R1, DSC$DayTwelve_R2 , DSC$DayTwelve_R3) , na.rm = TRUE)


Agenes <- DSCM[clustA,]
Bgenes <- DSCM[clustB,]
Cgenes <- DSCM[clustC,]
Dgenes <- DSCM[clustD,]
Egenes <- DSCM[clustE,]
Fgenes <- DSCM[clustF,]
Ggenes <- DSCM[clustG,]
Hgenes <- DSCM[clustH,]
Igenes <- DSCM[clustI,]

mygenes <- list(Agenes,Bgenes,Cgenes,Dgenes,Egenes,Fgenes,Ggenes,Hgenes, Igenes)


#==================================== plot function =========================================
plot.clusterRici<- function(w,x,y,z,xlab=NULL,title=title){
  w<- as.matrix(w)
  x<- as.matrix(x)
  y<- as.matrix(y)
  z<- as.matrix(z)
  sd.w<- apply(w,2,sd)
  sd.x<- apply(x,2,sd)
  sd.y<- apply(y,2,sd)
  sd.z<- apply(z,2,sd)
  mean.w<- apply(w,2,mean)
  mean.x<- apply(x,2,mean)
  mean.y<- apply(y,2,mean)
  mean.z<- apply(z,2,mean)
  
  length.x<- length(sd.w)+length(sd.x)+length(sd.y)+length(sd.z)
  
  # 	print(sd.x)
  # 	print(sd.y)
  # 	print(sd.z)
  #   print(mean.x)
  # 	print(mean.y)
  # 	print(mean.z)
  
  
  max.y<- ceiling(max(mean.w,mean.x,mean.y,mean.z) + max(sd.x,sd.y,sd.z))
  min.y<- floor(min(mean.w,mean.x,mean.y,mean.z) - max(sd.x,sd.y,sd.z))
  
  
  plot.new()
  plot.window(c(0,length.x+0.5 ),c(min.y,max.y))
  
  if(!is.null(xlab)){
    axis(1, at=seq(0, length(xlab)-1, 1),labels=xlab)
  }
  else{
    print("No x-axis labeling")
  }
  
  axis(2, at=seq(min.y, max.y, 1),cex.axis=1,las=2)
  #text(1,max.y,labels="rigid")
  #text(5,max.y,labels="critical")
  
  if(!is.null(title)){
    title(main=title)
  }
  
  points(seq(0.5, length.x-0.5, 1),c(mean.w,mean.x,mean.y,mean.z))
  lines(seq(0.5, length.x-0.5, 1),c(mean.w,mean.x,mean.y,mean.z))
  
  # 	lines(seq(0, length.x/3-1, 1),mean.x)
  # 	lines(seq(length.x/3,(length.x/3*2)-1, 1),mean.y)
  # 	lines(seq(length.x/3*2,length.x-1, 1),mean.z)
  
  
  plotbars(0.5,mean.w,sd.w)
  plotbars(1.5,mean.x,sd.x)
  plotbars(2.5,mean.y,sd.y)
  plotbars(3.5,mean.z,sd.z)
  
  #  	box()
}


plotbars<- function(x,y,sd,w=0.1){
  
  X<- c(x,x)
  Y<- c(y+sd,y-sd)
  lines(X,Y)
  X<- c(x-w,x+w)
  Y<- c(y+sd,y+sd)
  lines(X,Y)
  Y<- c(y-sd,y-sd)
  lines(X,Y)
}
 #==============================================================================================================================
pdf(paste(path,"LinesPlot_allScale2.pdf",sep = "/"), height=4,width=5)
par(mfrow=c(1,1))


for (i in (1:9)){
  mychar <- toupper(letters[i])
  mym2 <- mygenes[[i]]
  title = paste("Cluster",mychar)
  mym<- t(apply(as.matrix(mym2), 1, scale))
  w <- mym[,1]
  x <- mym[,2]
  y <- mym[,3]
  z <- mym[,4]
  
  plot.clusterRici(w,x,y,z,title = title)
}

dev.off()


boxplot(mym, ylim=c(0,9000))


