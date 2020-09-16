#---------------------------------------------------------------------
# FILE     : tnbctilsSurvIndivMETABRIC.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-15
# COMMENTS : This file reads in the survival and cibersort data and 
#            will generate cell specific histograms (cut offs found
#            from tnbctilsSurvMETABRIC.R)
#---------------------------------------------------------------------

library(ggplot2)

surFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcSurv.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
sur$OS_STATUS2 <- ifelse(sur$OS_STATUS == "0:LIVING", 0, 1)

vals <- "abs"
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\intensity\\"
dirSave2 <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\"

# 5 = CD8, 8 = CD4, 11 = gamma delta
cell_idx <- 5

cibersort <- rel[,c(1,cell_idx)]
cellType <- colnames(cibersort)[2]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)

sur <- sur[!(is.na(sur$sum)),]

#---------
# ggplot specific histo

# CD8
low <- 0.0613 #0.1011
high <- 0.1023 #0.1726
legtitle <- "CD8 T Cells"
axis <- "Quantity of CD8 T Cells"
pos <- c(0.7, 0.75)

# median
low <- 0.05809401 #0.0849
high <- 0.05809401  #0.0849
legtitle <- "CD8 T Cells"
axis <- "Quantity of CD8 T Cells"
pos <- c(0.7, 0.75)

# T Cells gamma delta
low <- 0.0382
high <- 0.0662
legtitle <- "T Cells gamma delta"
axis <- "Quantity of T Cells gamma delta"
pos <- c(0.7, 0.75)

# median
low <- 0.03618013
high <- 0.03618013
legtitle <- "T Cells gamma delta"
axis <- "Quantity of T Cells gamma delta"
pos <- c(0.7, 0.75)

# CD4
#low <- 0.0142 #0.0348
#high <- 0.0325 #0.0646
#legtitle <- "CD4 Memory\nActivated T Cells"
#axis <- "Quantity of CD4\nMemory Activated T Cells"
#pos <- c(0.7, 0.75)

# median
#low <- 0.007955437
#high <- 0.007955437
#legtitle <- "CD4 Memory\nActivated T Cells"
#axis <- "Quantity of CD4\nMemory Activated T Cells"
#pos <- c(0.7, 0.75)

#-----

h<-hist(cibersort[,2], breaks=10, plot=F)
cibersort$range <- ifelse(cibersort[,2] < low, "low", ifelse(cibersort[,2] > high, "high", "medium"))
cibersort$range <- factor(cibersort$range, levels=c("low", "medium", "high"), ordered=T)

# median
h<-hist(cibersort[,2], breaks=10, plot=F)
cibersort$range <- ifelse(cibersort[,2] < low, "low", ifelse(cibersort[,2] > high, "high", "medium"))
cibersort$range <- factor(cibersort$range, levels=c("low", "high"), ordered=T)

# plot using density
#  ggplot(, aes(x=cibersort[,2])) +
#    geom_histogram(aes(y=..density.., fill=cibersort$range), breaks=h$breaks, colour="black") +
#    scale_fill_manual(values = c("#F8766D", "grey", "#00BFC4")) +
#    geom_vline(aes(xintercept = 0.1011), colour="black", size=1) +
#    geom_vline(aes(xintercept = 0.1726), colour="black", size=1) +
#    stat_function(fun = dnorm, args = c(mean = mean(cibersort[,2]), sd = sd(cibersort[,2])))

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "hist.", vals, ".", cellType, ".2", ".png"), 
    width=width, height=height, res=300)

ggplot(cibersort,aes(x=cibersort[,2])) +
  geom_histogram(aes(fill=range), breaks=h$breaks, colour="black") +
  scale_fill_manual(values = c("#F8766D", "grey", "#00BFC4")) +
  geom_vline(aes(xintercept = low), colour="black", size=1) +
  geom_vline(aes(xintercept = high), colour="black", size=1) +
  stat_function(fun = function(x, mean, sd, n){ 0.05*n*dnorm(x=x, mean=mean, sd=sd)}, 
                args=c(mean=mean(cibersort[,2]), sd=sd(cibersort[,2]), n=length(cibersort[,2])),
                size=1) +
  #ggtitle("Histogram of CD8 T Cells") +
  xlab(axis) +
  ylab("Count of Samples") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(legend.text = element_text(size=20, colour="black")) +  
  theme(legend.title = element_text(size=22, colour="black", face="bold")) +
  guides(fill=guide_legend(title=legtitle)) +
  theme(legend.position = pos)

dev.off()

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "hist.", vals, ".", cellType, ".3", ".png"), 
    width=width, height=height, res=300)

# median
ggplot(cibersort, aes(x=cibersort[,2])) +
  geom_histogram(aes(fill=range), breaks=h$breaks, colour="black") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  geom_vline(aes(xintercept = low), colour="black", size=1) +
  geom_vline(aes(xintercept = high), colour="black", size=1) +
  stat_function(fun = function(x, mean, sd, n){ 0.05*n*dnorm(x=x, mean=mean, sd=sd)}, 
                args=c(mean=mean(cibersort[,2]), sd=sd(cibersort[,2]), n=length(cibersort[,2])),
                size=1) +
  #ggtitle("Histogram of CD8 T Cells") +
  xlab(axis) +
  ylab("Count of Samples") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(legend.text = element_text(size=20, colour="black")) +  
  theme(legend.title = element_text(size=22, colour="black", face="bold")) +
  guides(fill=guide_legend(title=legtitle)) +
  theme(legend.position = pos)

dev.off()

# ggplot specific histo
#---------

# temporary code for saving data
save <- sur
save$imm <- NA
save[save$sum < low & !is.na(save$sum),]$imm <- "low"
save[save$sum > high & !is.na(save$sum),]$imm <- "high"
toFile <- data.frame(barcode=save$PATIENT_ID, imm=save$imm)
write.table(toFile, paste0(dirSave2, "os.abs.T.cells.CD8.0.25sd.txt"), sep="\t", row.names=F, quote=F)
write.table(toFile, paste0(dirSave2, "os.abs.T.cells.CD8.mediansd.txt"), sep="\t", row.names=F, quote=F)
write.table(toFile, paste0(dirSave2, "os.abs.T.cells.gamma.delta.0.25sd.txt"), sep="\t", row.names=F, quote=F)
write.table(toFile, paste0(dirSave2, "os.abs.T.cells.gamma.delta.mediansd.txt"), sep="\t", row.names=F, quote=F)
write.table(toFile, paste0(dirSave2, "os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), sep="\t", row.names=F, quote=F)
write.table(toFile, paste0(dirSave2, "os.abs.T.cells.CD4.memory.activated.mediansd.txt"), sep="\t", row.names=F, quote=F)

write.table(sur, paste0(dirSave2, "os.abs.T.cells.CD8.0.25sd.survival.txt"), sep="\t", row.names=F, quote=F)
write.table(sur, paste0(dirSave2, "os.abs.T.cells.gamma.delta.0.25sd.survival.txt"), sep="\t", row.names=F, quote=F)
write.table(sur, paste0(dirSave2, "dfs.abs.T.cells.CD4.memory.activated.0.25sd.survival.txt"), sep="\t", row.names=F, quote=F)
