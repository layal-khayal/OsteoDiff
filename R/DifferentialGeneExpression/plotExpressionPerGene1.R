# this script to plot the expression of candidate genes

library(ggplot2)
library(grid)
library(gridExtra)
library(DESeq2)
library(scales)

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


#=================== plotting functions =======================================
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  theme(plot.title = element_text(face = "bold",size = 14, hjust = 0, family = 'sans'),
        text = element_text(size = 11),
        panel.background = element_blank(),
        axis.title = element_text(face = "bold",size = 13, family = 'sans'),
        axis.text.x = element_text(size = 12, colour = 'black', family = 'sans'),
        axis.text.y = element_text(size = 12, colour = 'black', family = 'sans'),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        legend.key = element_rect(colour = 'white'),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(face= "italic" ,size = 12, colour = 'black', family = 'sans')
        #plot.margin = margin(0.2,1,1,0.2,"cm")
  )
  
  
  
  
          #legend.key.size= unit(0.2, "cm"),
          #legend.title = element_text(face="italic"),
          #legend.title=element_blank(),
          #plot.margin=unit(c(10,5,5,5),"mm"),
         
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#=================================================================================
#-------------------Bmp2
Bmp2<- plotCounts(ddseq, gene= "Bmp2", intgroup = "condition", normalized = TRUE, transform = FALSE, returnData = TRUE)
BmpMean <- plyr::ddply(Bmp2, "condition", plyr::summarise, mean = mean(count))
g1<-ggplot(Bmp2, aes(x=condition, y=count))+ 
  geom_point(aes(shape="Normalized counts", colour= "Replicates"), size=3, position = position_jitter(height = 0.01, width = 0.1) )+
  geom_point(data = BmpMean, aes(x= condition, y = mean, shape="Mean count",colour="Mean of replicates"), size = 4)+
  ggtitle("Bmp2")+
  scale_shape_manual(values=c(4, 17),
                     breaks=c("Normalized counts", "Mean count" ),
                     name="") +
  scale_color_manual (values = c("red","blue"),
                      breaks = c("Replicates", "Mean of replicates" ),
                      name="") +
  
  scale_y_continuous(limits=c(0,max(Bmp2$count)*1.1),
                     labels = scientific_format(digits = 2) ) +
  theme_Publication()


#--------------- Runx2
Runx<- plotCounts(ddseq, gene= "Runx2", intgroup = "condition", normalized = TRUE, transform = FALSE, returnData = TRUE)
RunMean <- plyr::ddply(Runx, "condition", plyr::summarise, mean = mean(count))
g2<-ggplot(Runx, aes(x=condition, y=count))+ 
  geom_point(shape=17,colour="blue", size=3, position = position_jitter(height = 0.01, width = 0.1) )+
  geom_point(data = RunMean, aes(x= condition, y = mean),shape= 4 ,colour="red", size = 4)+
  ggtitle("Runx2")+
  scale_y_continuous(limits=c(0,max(Runx$count)*1.1),
                     labels = scientific_format(digits = 2) ) +
  theme_Publication()

# --------------- col1a1
col1<- plotCounts(ddseq, gene= "Col1a1", intgroup = "condition", normalized = TRUE, transform = FALSE, returnData = TRUE)
col1Mean <- plyr::ddply(col1, "condition", plyr::summarise, mean = mean(count))
g3<-ggplot(col1, aes(x=condition, y=count))+ 
  geom_point(shape=17,colour="blue", size=3, position = position_jitter(height = 0.01, width = 0.1) )+
  geom_point(data = col1Mean, aes(x= condition, y = mean),shape= 4 ,colour="red", size = 4)+
  ggtitle("Col1a1")+
  scale_y_continuous(limits=c(0,max(col1$count)*1.1),
                     labels = scientific_format(digits = 2) ) +
  theme_Publication()

#----------------Mmp13
Mmp<- plotCounts(ddseq, gene= "Mmp13", intgroup = "condition", normalized = TRUE, transform = FALSE, returnData = TRUE)
MmpMean <- plyr::ddply(Mmp, "condition", plyr::summarise, mean = mean(count))

g4<-ggplot(Mmp, aes(x=condition, y=count))+ 
  geom_point(shape=17,colour="blue", size=3, position = position_jitter(height = 0.01, width = 0.1) )+
  geom_point(data = MmpMean, aes(x= condition, y = mean),shape= 4 ,colour="red", size = 4)+
  ggtitle("Mmp13")+
  scale_y_continuous(limits=c(0,max(Mmp$count)*1.1),
                     labels = scientific_format(digits = 2) ) +
  theme_Publication()

# ------------ Runx2
runx <- plotCounts(ddseq, gene= "Runx2", intgroup = "condition", normalized = TRUE, transform = FALSE, returnData = TRUE)
runxMean <- plyr::ddply(runx, "condition", plyr::summarise, mean = mean(count))
g5 <- ggplot(runx, aes(x=condition, y=count))+ 
  geom_point(shape=17,colour="blue", size=3, position = position_jitter(height = 0.01, width = 0.1) )+
  geom_point(data = runxMean, aes(x= condition, y = mean),shape= 4 ,colour="red", size = 4)+
  ggtitle("Runx2")+
  theme_Publication()

#--------------------- Sost
Sost<- plotCounts(ddseq, gene= "Sost", intgroup = "condition", normalized = TRUE, transform = FALSE, returnData = TRUE)
SostMean <- plyr::ddply(Sost, "condition", plyr::summarise, mean = mean(count))
g6<-ggplot(Sost, aes(x=condition, y=count))+ 
  geom_point(shape=17,colour="blue", size=3, position = position_jitter(height = 0.01, width = 0.1) )+
  geom_point(data = SostMean, aes(x= condition, y = mean),shape= 4 ,colour="red", size = 4)+
  ggtitle("Sost")+
  theme_Publication()



#====================================================== plot common legend ================================================
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right"),
                                       labs=list(), labpos=list(c(0.5,0.02), c(0.02,0.5))) {
  
  plots <- list(...)
  # for position
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none", axis.title = element_blank()))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight),
                                            vp=viewport(width=0.9, height=0.9)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  if(!length(labs) == 0){
    grid.text(labs[1], x=labpos[[1]][1], y=labpos[[1]][2], gp=gpar(fontsize=16, fontfamily= 'Times', fontface= c('bold','italic')))
    grid.text(labs[2], x=labpos[[2]][1], y=labpos[[2]][2], rot=90, gp=gpar(fontsize=16, fontfamily= 'Times', fontface= c('bold','italic')))
  }

  
  # return gtable invisibly
  invisible(combined)
  
}

grid_arrange_shared_legend(g1, g2, g3, g4,g5,g6, ncol= 2,nrow = 3, labs=list("Differentiation Time-point", "Expression / normalized counts"))
#===============================================================================================================================

pdf("GE_counts.pdf", width = 6.68, height = 7.2, pointsize = 16, family = 'Times', onefile=FALSE)
grid_arrange_shared_legend(g1, g2, g3, g4,g5,g6, ncol= 2,nrow = 3, labs=list("Differentiation Time-point", "Expression / normalized counts"))
dev.off()

  
     
        
labpos=list(c(0.5,0.03), c(0.03,0.5))
grid.arrange(g1,g2, g3, g4, ncol=2)
grid.text("Lolo", x=0.5, y=0.03, gp=gpar(fontsize=16))
grid.text("Shatora", x=0.03, y=0.5, rot=90, gp=gpar(fontsize=16))

