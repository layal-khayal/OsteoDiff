install.packages("ggplot2")
install.packages("calibrate")
library("ggplot2")
library("reshape2")

rootout = "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/Quality_control/Python/coverage/"

#in Ubunto my computer

rootinput= "/media/layal/Ubuntu_data/Osteo_project/Osteo_2017/Quality_control/Python/coverage/coverage_output/" 

OsteoSample <- c("Day_0_R1", "Day_0_R2", "Day_0_R3", 
                 "Day_3_R1", "Day_3_R2", "Day_3_R3",
                 "Day_6_R1", "Day_6_R2", "Day_6_R3",
                 "Day_12_R1", "Day_12_R2", "Day_12_R3")

Samp_Cov <- c("_exon_cov_66604.txt", "_exon_cov_66605.txt", "_exon_cov_66606.txt",
              "_exon_cov_66607.txt", "_exon_cov_66608.txt", "_exon_cov_66609.txt",
              "_exon_cov_66610.txt", "_exon_cov_66611.txt", "_exon_cov_66613.txt",
              "_exon_cov_66614.txt", "_exon_cov_66615.txt", "_exon_cov_66616.txt") 

sam_key <- c("first", "middle","last")

cond <- c("0", "3", "6", "12")


for(k in sam_key){
  j=0
  
  X= paste("DaysAvg",k,sep = "_")
  
  #iterate over 12 samples, and every 3 samples is a condition j
  for(i in seq(1,12,3)){
    
    Rep_1 <- assign(paste( OsteoSample[i], k,sep = "_"), read.table(paste(rootinput, k , Samp_Cov[i], sep = "")))
    
     
    Rep_2 <-assign(paste( OsteoSample[i+1], k,sep = "_"), read.table(paste(rootinput, k , Samp_Cov[i+1], sep = "")))
    
     
    Rep_3 <- assign(paste( OsteoSample[i+2], k,sep = "_"), read.table(paste(rootinput, k , Samp_Cov[i+2], sep = "")))
    
    j=j+1
    
    cov_con <- assign(paste("Cov",cond[j],k,sep = "_"), cbind(Rep_1[[2]], Rep_2[[2]] , Rep_3[[2]]))
    
     
    avg_con <- assign(paste("Avg",cond[j], k, sep = "_"), rowMeans(cov_con))
    colNam<- paste("Day",cond[j],sep = "_")
    
    if (j==1){
      df <- data.frame(avg_con)
      colnames(df) <- colNam
    }
    else{
    df[[colNam]]<-avg_con 
    }
    assign(X, df)
  }
}

#  ----------------====================== Plotting ==========================----------------------------------

# ------------------ plot coverage first exon
#we would like is to have our data in 2 columns, where the first column contains the data values, 
#and the second column contains the sample name

dataDenFst<- stack(DaysAvg_first)

head(dataDenFst)
head(log(dataDenFst$values))

dataFst <- dataDenFst[which(dataDenFst$values > 0),]
colnames(dataFst)<-c("values","Conditions")
head(dataFst)
head(log(dataFst$values))

options(digits=9)
pdf( file=paste(rootout,"firstExonsCovDensity.pdf",sep = ""), width=12, height=10, pointsize=12)

ggplot(dataFst , aes(x =values)) + 
  geom_density(aes(group=Conditions, colour=Conditions, fill=Conditions), alpha=0.5) +
  scale_x_log10(limits=c(0.1, 10000))+ 
  xlab("Mean of base coverage")+
  ggtitle("Coverage Density of First Exons")+ 
  theme(plot.title = element_text(hjust = 0.5),text= element_text(size = 25))
dev.off()


# ------------------ plot coverage second exon ---------------------------------------------------------------------------

dataDenMid<- stack(DaysAvg_middle)

dataMid <- dataDenMid[which(dataDenMid$values > 0),]
colnames(dataMid)<-c("values","Conditions")

pdf( file=paste(rootout,"MiddleExonsCovDensity.pdf",sep = ""), width=12, height=10, pointsize=12)

ggplot(dataMid , aes(x =values)) + 
  geom_density(aes(group=Conditions, colour=Conditions, fill=Conditions), alpha=0.5) +
  scale_x_log10(limits=c(0.1, 10000))+ 
  xlab("Mean of base coverage")+
  ggtitle("Coverage Density of Middle Exons")+ 
  theme(plot.title = element_text(hjust = 0.5),text= element_text(size = 25))
dev.off()

# ------------------ plot coverage last exon ---------------------------------------------------------------------------

dataDenlas<- stack(DaysAvg_last)

datalas <- dataDenlas[which(dataDenlas$values > 0),]
colnames(datalas)<-c("values","Conditions")

pdf( file=paste(rootout,"LastExonsCovDensity.pdf",sep = ""), width=12, height=10, pointsize=12)

ggplot(datalas , aes(x =values)) + 
  geom_density(aes(group=Conditions, colour=Conditions, fill=Conditions), alpha=0.5) +
  scale_x_log10(limits=c(0.1, 10000))+ 
  xlab("Mean of base coverage")+
  ggtitle("Coverage Density of Last Exons")+ 
  theme(plot.title = element_text(hjust = 0.5),text= element_text(size = 25))
dev.off()





















pdf( file=paste(rootout,"firstExonsCovDensity_1.pdf",sep = ""), width=12, height=10, pointsize=12)
ggplot(dataFst , aes(x =values)) + geom_density(aes(group=Conditions, colour=Conditions, fill=Conditions), alpha=0.5) +scale_x_log10()+ xlab("Mean of base coverage")+
  ggtitle("Coverage Density of First Exons")+ theme(plot.title = element_text(hjust = 0.5), text= element_text(size = 20))
dev.off()


