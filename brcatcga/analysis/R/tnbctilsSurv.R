#---------------------------------------------------------------------
# FILE     : tnbctilsSurv.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-08-05
# COMMENTS : read in a tab-delimited file with TCGA info, including 
#            survival information, and summarize it for analysis by R
#---------------------------------------------------------------------

library(ggplot2)
library(survival)
library(GGally)
library(ggfortify)

# cibersort
#rel <- read.delim("~/brca/analysis/ANALYSIS/CIBERSORT.Output_Job4.txt", header=T, stringsAsFactors=F)
#rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

##surFile <- "~/brca/brcatcga/analysis/R/tnbc.tils.append.surv.idc.txt"
surFile <- "F:\\TNBC TILS\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"

vals <- "abs"
##rel <- read.delim("~/brca/brcatcga/analysis/R/cibersort/data/CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

#for (cell_idx in 2:(dim(rel)[2]-3)){
#for (cell_idx in 3:3){
# 5 = CD8, 8 = CD4
  cell_idx <- 5
  cibersort <- rel[,c(1,cell_idx)]
  cellType <- colnames(cibersort)[2]
  colnames(cibersort)[2] <- "sum"

  sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
  #cibersort <- read.delim("~/brca/analysis/R/tnbcGeneExpCiberSort.htseq.txt", header=T, stringsAsFactors=F, na.strings=(""))
  sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
  #smpl.surv <- Surv(sur$os_days, sur$os_status)~1
  sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
  sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
  sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
  sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
  #sur <- sur[!(is.na(sur$pathologic_stage)),]
  sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)

  # analysis of cibersort
  sur <- sur[!(is.na(sur$sum)),]

  x <- sur$sum

  width <- 600
  #ratio <- 461/572
  ratio <- 547/600
  height <- trunc(width*ratio)
  png(paste0("~/brca/brcatcga/analysis/R/cibersort/fpkm/hist.", vals, ".", cellType, ".png"), width=width, height=height)
  
  # automated histo
  h<-hist(x, breaks=10, col="red", xlab="quantity of cell type",
          main=paste0("Histogram of ", cellType))
  xfit<-seq(min(x),max(x),length=40)
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
  yfit <- yfit*diff(h$mids[1:2])*length(x)
  lines(xfit, yfit, col="blue", lwd=2)
  
  dev.off()
  
  #---------
  # ggplot specific histo

  # CD8
  low <- 0.1011
  high <- 0.1726
  legtitle <- "CD8 T Cells"
  axis <- "Quantity of CD8 T Cells"
  pos <- c(0.7, 0.75)
  
  # median
  low <- 0.0849
  high <- 0.0849
  legtitle <- "CD8 T Cells"
  axis <- "Quantity of CD8 T Cells"
  pos <- c(0.7, 0.75)
  
  # CD4
  low <- 0.0348
  high <- 0.0646
  legtitle <- "CD4 Memory\nActivated T Cells"
  axis <- "Quantity of CD4\nMemory Activated T Cells"
  pos <- c(0.7, 0.75)
  
  # median
  low <- 0.0258
  high <- 0.0258
  legtitle <- "CD4 Memory\nActivated T Cells"
  axis <- "Quantity of CD4\nMemory Activated T Cells"
  pos <- c(0.7, 0.75)
  
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
  #png(paste0("~/brca/brcatcga/analysis/R/cibersort/fpkm/hist.", vals, ".", cellType, ".2", ".png"), width=width, height=height)
  png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\cibersort\\fpkm\\hist.", vals, ".", cellType, ".2", ".png"), 
      width=width, height=height, res=300)
  
  ggplot(, aes(x=cibersort[,2])) +
    geom_histogram(aes(fill=cibersort$range), breaks=h$breaks, colour="black") +
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
  #png(paste0("~/brca/brcatcga/analysis/R/cibersort/fpkm/hist.", vals, ".", cellType, ".2", ".png"), width=width, height=height)
  png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\cibersort\\fpkm\\hist.", vals, ".", cellType, ".3", ".png"), 
      width=width, height=height, res=300)
  
  # median
  ggplot(, aes(x=cibersort[,2])) +
    geom_histogram(aes(fill=cibersort$range), breaks=h$breaks, colour="black") +
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

  m <- mean(sur$sum, na.rm=T)
  sd <- sd(sur$sum, na.rm=T)
  #for (times in c(0.25, 0.5, 1)){
  times <- 0.25
  lowEnd <- round(m - (times * sd), digits=4)
  #if (lowEnd <= 0){ next }
  highEnd <- round(m + (times * sd), digits=4)
  
  # median
  m <- mean(sur$sum, na.rm=T)
  med <- median(sur$sum, na.rm=T)
  m <- med
  sd <- sd(sur$sum, na.rm=T)
  #for (times in c(0.25, 0.5, 1)){
  #times <- 0.25
  times <- "median"
  lowEnd <- m
  highEnd <- m

    # low and high
    #lowEnd <- 20
    #highEnd <- 40

  # temporary code for saving data
  #save <- sur
  #save$imm <- NA
  #save[save$sum < lowEnd & !is.na(save$sum),]$imm <- "low"
  #save[save$sum > highEnd & !is.na(save$sum),]$imm <- "high"
  #toFile <- data.frame(barcode=save$bcr_patient_barcode, imm=save$imm)
  #write.table(toFile, "~/brca/brcatcga/analysis/R/os.abs.T.cells.CD8.0.25sd.txt", sep="\t", row.names=F, quote=F)
  #write.table(toFile, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.mediansd.txt", sep="\t", row.names=F, quote=F)
  #write.table(toFile, "~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.txt", sep="\t", row.names=F, quote=F)
  #write.table(toFile, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.mediansd.txt", sep="\t", row.names=F, quote=F)
  
    # lowEnd, highEnd
    sur$imm <- NA
    sur[sur$sum < lowEnd & !is.na(sur$sum),]$imm <- paste0("<", toString(lowEnd))
    sur[sur$sum > highEnd & !is.na(sur$sum),]$imm <- paste0(">", toString(highEnd))
    sur <- sur[!(is.na(sur$imm)),]
    sur$imm <- factor(sur$imm, levels=c(paste0("<", toString(lowEnd)), paste0(">", toString(highEnd))), ordered=T)
    lowNum <- dim(sur[sur$imm==paste0("<", toString(lowEnd)),])[1]
    highNum <- dim(sur[sur$imm==paste0(">", toString(highEnd)),])[1]
    
    # mean/median
    m <- med
    sur$imm <- NA
    sur[sur$sum < m & !is.na(sur$sum),]$imm <- paste0("<", toString(m))
    sur[sur$sum > m & !is.na(sur$sum),]$imm <- paste0(">", toString(m))
    sur <- sur[!(is.na(sur$imm)),]
    sur$imm <- factor(sur$imm, levels=c(paste0("<", toString(m)), paste0(">", toString(m))), ordered=T)
    lowNum <- dim(sur[sur$imm==paste0("<", toString(m)),])[1]
    highNum <- dim(sur[sur$imm==paste0(">", toString(m)),])[1]

    #write.table(sur, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.survival.txt", sep="\t", row.names=F, quote=F)
    #write.table(sur, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\dfs.abs.T.cells.CD4.memory.activated.0.25sd.survival.txt", sep="\t", row.names=F, quote=F)
    
    # km
    smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$imm
    sigTest <- survdiff(smpl.surv)
    pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
    smpl.fit <- survfit(smpl.surv)

    width <- 600
    #ratio <- 461/572
    ratio <- 547/600
    height <- trunc(width*ratio)
    png(paste0("~/brca/brcatcga/analysis/R/cibersort/fpkm/os.", vals, ".", cellType, ".", times, "sd.png"), width=width, height=height)

    palette <- c("#F8766D", "#00BFC4")
    plot(smpl.fit, lty=c(1,1), main=paste0("Disease Free Survival ", cellType, "\nPvalue =", round(pvalue, digits=3)),
         xlab="days", ylab="survival", col=palette, lwd=2, mark=108, cex=1.5)
    legend(500,0.3, legend=c(paste0(">", toString(highEnd), " ", cellType, " n=", highNum), paste0("<", toString(lowEnd), " ", cellType, " n=", lowNum)), 
           lty=c(1,1), col=rev(palette), lwd=2)

    dev.off()
  
    #---
    # fancy plot
  # os
  width <- 2500
  #ratio <- 461/572
  ratio <- 547/600
  height <- trunc(width*ratio)
  #png(paste0("~/brca/brcatcga/analysis/R/cibersort/fpkm/os.", vals, ".", cellType, ".", times, "sd.2.png"), width=width, height=height)
  png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\cibersort\\fpkm\\os.", vals, ".", cellType, ".", times, "sd.2.png"), 
      width=width, height=height, res=300)
  
    palette <- c("#F8766D", "#00BFC4")
    par(mar=c(5.1,5.1,4.1,2.1))
    plot(smpl.fit, lty=c(1,1),
         xlab="Days", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
    legend(500,0.35, legend=c(paste0("High, n=", highNum), paste0("Low, n=", lowNum)), 
           lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)
    
    dev.off()
  
  # dfs
  width <- 2500
  #ratio <- 461/572
  ratio <- 547/600
  height <- trunc(width*ratio)
  #png(paste0("~/brca/brcatcga/analysis/R/cibersort/fpkm/dfs.", vals, ".", cellType, ".", times, "sd.2.png"), width=width, height=height)
  png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\cibersort\\fpkm\\dfs.", vals, ".", cellType, ".", times, "sd.2.png"), 
      width=width, height=height, res=300)
  
  palette <- c("#F8766D", "#00BFC4")
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(smpl.fit, lty=c(1,1),
       xlab="Days", ylab="Disease Free Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  legend(500,0.35, legend=c(paste0("High, n=", highNum), paste0("Low, n=", lowNum)), 
         lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)
  
  dev.off()
  
  # fancy plot
  #---
 
  #---------------------
  # analysis of continuous
  #----------------------
  # cox
  cox <- coxph(Surv(os_days, os_status)~sum, data=sur)
  summary(cox)
  
  cox <- coxph(Surv(os_days, os_status)~pathologic_stage, data=sur)
  summary(cox)
  
  cox <- coxph(Surv(os_days, os_status)~age_at_initial_pathologic_diagnosis, data=sur)
  summary(cox)
  
  cox <- coxph(Surv(os_days, os_status)~sum+pathologic_stage, data=sur)
  summary(cox)
  
  # cox
  cox <- coxph(Surv(dfs_days, dfs_status)~sum, data=sur)
  summary(cox)
  
  cox <- coxph(Surv(dfs_days, dfs_status)~pathologic_stage, data=sur)
  summary(cox)
  
  cox <- coxph(Surv(dfs_days, dfs_status)~age_at_initial_pathologic_diagnosis, data=sur)
  summary(cox)
  
  cox <- coxph(Surv(dfs_days, dfs_status)~sum+pathologic_stage, data=sur)
  summary(cox)
  
#  }
#}

sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))

#smpl.surv <- Surv(sur$os_days, sur$os_status)~1

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III"
#sur <- sur[!(is.na(sur$pathologic_stage)),]
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I", "Stage II", "Stage III", "Stage IV"))


ggsurv(smpl.fit, xlab="days", size.est=0.75, cens.size=4, cens.shape=3)

# cox
cox <- coxph(Surv(os_days, os_status)~imm, data=sur)
summary(cox)
df <- with(sur, data.frame(tils = c("<20%", ">20%")))

cox <- coxph(Surv(os_days, os_status)~pathologic_stage, data=sur)
summary(cox)
df <- with(sur, data.frame(pathologic_stage = c("Stage I", "Stage II", "Stage III", "Stage IV")))

cox <- coxph(Surv(os_days, os_status)~imm+pathologic_stage, data=sur)
summary(cox)
df <- with(sur, data.frame(tils = c("<20%", ">20%"), pathologic_stage = rep(levels(pathologic_stage)[2], 2)))

cox <- coxph(Surv(os_days, os_status)~imm+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~imm+age_at_initial_pathologic_diagnosis+pathologic_stage, data=sur)
summary(cox)

#---------------------
# analysis of tils
#----------------------
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
sur <- sur[!(is.na(sur$tils)),]
#sur <- sur[!(is.na(sur$pathologic_stage)),]
#sur <- sur[sur$pathologic_stage=="Stage II",]
sur$tils <- factor(sur$tils, levels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"), ordered=T)

sur[sur$tils %in% c("<1%"),]$tils <- 0
sur[sur$tils %in% c("1-10%"),]$tils <- 1
sur[sur$tils %in% c("10-20%"),]$tils <- 2
sur[sur$tils %in% c("20-30%"),]$tils <- 3
sur[sur$tils %in% c("30-40%"),]$tils <- 4
sur[sur$tils %in% c("40-50%"),]$tils <- 5
sur[sur$tils %in% c("50-60%"),]$tils <- 6
sur[sur$tils %in% c(">70%"),]$tils <- 7
sur$tils <- as.numeric(sur$tils)

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
#sur[sur$pathologic_stage %in% c("Stage I-II"),]$pathologic_stage <- 0
#sur[sur$pathologic_stage %in% c("Stage III-IV"),]$pathologic_stage <- 1
#sur[sur$pathologic_stage %in% c("Stage III-IV"),]$pathologic_stage <- 2
#sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- 2
#sur$pathologic_stage <- as.numeric(sur$pathologic_stage)

#sur[sur$age_at_initial_pathologic_diagnosis <= 50,]$age_at_initial_pathologic_diagnosis <- 0
#sur[sur$age_at_initial_pathologic_diagnosis >50,]$age_at_initial_pathologic_diagnosis <- 1

# <1%, 1-20%, 20-40%, >40%
#sur[sur$tils %in% c("<1%"),]$tils <- "<1%"
#sur[sur$tils %in% c("1-10%", "10-20%"),]$tils <- "1-20%"
#sur[sur$tils %in% c("20-30%", "30-40%"),]$tils <- "20-40%"
#sur[sur$tils %in% c("40-50%", "50-60%", "60-70%", ">70%"),]$tils <- ">40%"
#sur$tils <- factor(sur$tils, levels=c("<1%", "1-20%", "20-40%", ">40%"))

# <20% and >20%
sur[sur$tils %in% c("<1%", "1-10%", "10-20%"),]$tils <- "<20%"
sur[sur$tils %in% c("20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"),]$tils <- ">20%"
sur$tils <- factor(sur$tils, levels=c("<20%", ">20%"))
lowNum <- dim(sur[sur$tils=="<20%",])[1]
highNum <- dim(sur[sur$tils==">20%",])[1]

# <30% and >30%
sur[sur$tils %in% c("<1%", "1-10%", "10-20%", "20-30%"),]$tils <- "<30%"
sur[sur$tils %in% c("30-40%", "40-50%", "50-60%", "60-70%", ">70%"),]$tils <- ">30%"
sur$tils <- factor(sur$tils, levels=c("<30%", ">30%"), ordered=T)
lowNum <- dim(sur[sur$tils=="<30%",])[1]
highNum <- dim(sur[sur$tils==">30%",])[1]

#sur[sur$tils %in% c("<1%", "1-10%", "10-20%"),]$tils <- "<20%"
#sur[sur$tils %in% c("30-40%", "40-50%", "50-60%", "60-70%", ">70%"),]$tils <- ">30%"
#sur$tils <- factor(sur$tils, levels=c("<20%", ">30%"))

#sur$tils <- factor(sur$tils, levels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", ">70%"))

# cox
cox <- coxph(Surv(os_days, os_status)~tils, data=sur)
summary(cox)
df <- with(sur, data.frame(tils = c("<30%", ">30%")))

cox <- coxph(Surv(os_days, os_status)~pathologic_stage, data=sur)
summary(cox)
df <- with(sur, data.frame(pathologic_stage = c("Stage I", "Stage II", "Stage III", "Stage IV")))

cox <- coxph(Surv(os_days, os_status)~age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~tils+pathologic_stage, data=sur)
df <- with(sur, data.frame(tils = c("<30%", ">30%"), pathologic_stage = rep(levels(pathologic_stage)[2], 2)))
summary(cox)

cox <- coxph(Surv(os_days, os_status)~tils+age_at_initial_pathologic_diagnosis+pathologic_stage, data=sur)
cox <- coxph(Surv(os_days, os_status)~tils+pathologic_stage, data=sur)
cox <- coxph(Surv(os_days, os_status)~tils+strata(pathologic_stage), data=sur)
plot(survfit(cox))
summary(cox)
#df <- with(sur, data.frame(tils = c("<1%", "1-20%", "20-40%", ">40%"), pathologic_stage = rep(levels(pathologic_stage)[2], 2)))
df <- with(sur, data.frame(tils = c("<20%", ">20%"), age_at_initial_pathologic_diagnosis = rep(mean(age_at_initial_pathologic_diagnosis), 2), 
                                                                                               pathologic_stage = rep(levels(pathologic_stage)[2], 2)))

summary(survfit(cox, newdata=df))
plot(survfit(cox, newdata=df))

width <- 600
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png("~/brca/analysis/R/os_coxph_tils.stage_sII.png", width=width, height=height)

palette <- c("#F8766D", "#00BFC4")
plot(survfit(cox, newdata=df), lty=c(1,1), main=paste("Overall Survival"),
     xlab="days", ylab="survival", col=palette, lwd=2, mark=108, cex=1.5)
legend(500,0.3, legend=c(paste0(">20% TILs"), paste0("<20% TILs")), 
       lty=c(1,1), col=rev(palette), lwd=2)

dev.off()

# km os
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$tils
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
##png("~/brca/brcatcga/analysis/R/os2.png", width=width, height=height)
png("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os2-30.png", width=width, height=height, res=300)

#palette <- c("#F8766D", "#00BFC4")
#plot(smpl.fit, lty=c(1,1), main=paste("Overall Survival\nPvalue =", round(pvalue, digits=2)),
#     xlab="days", ylab="survival", col=palette, lwd=2, mark=108, cex=1.5)
#legend(500,0.3, legend=c(paste0(">20% TILs n=", highNum), paste0("<20% TILs n=", lowNum)), 
#       lty=c(1,1), col=rev(palette), lwd=2)

palette <- c("#F8766D", "#00BFC4")
par(mar=c(5.1,5.1,4.1,2.1))
plot(smpl.fit, lty=c(1,1),
     xlab="Days", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(500,0.35, legend=c(paste0(">30% TILs, n=", highNum), paste0("<30% TILs, n=", lowNum)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()

# km dfs
smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$tils
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
##png("~/brca/analysis/R/dfs2.png", width=width, height=height)
png("F:\\TNBC TILS\\brcatcga\\analysis\\R\\dfs2-30.png", 
    width=width, height=height, res=300)

#palette <- c("#F8766D", "#00BFC4")
#plot(smpl.fit, lty=c(1,1), main=paste("Disease Free Survival\nPvalue =", round(pvalue, digits=2)),
#     xlab="days", ylab="survival", col=palette, lwd=2, mark=108, cex=1.5)
#legend(500,0.3, legend=c(paste0(">20% TILs n=", highNum), paste0("<20% TILs n=", lowNum)), 
#       lty=c(1,1), col=rev(palette), lwd=2)

palette <- c("#F8766D", "#00BFC4")
par(mar=c(5.1,5.1,4.1,2.1))
plot(smpl.fit, lty=c(1,1),
     xlab="Days", ylab="Disease Free Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(500,0.35, legend=c(paste0(">30% TILs, n=", highNum), paste0("<30% TILs, n=", lowNum)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~tils, data=sur)
summary(cox)
df <- with(sur, data.frame(tils = c("<20%", ">20%")))

cox <- coxph(Surv(dfs_days, dfs_status)~pathologic_stage, data=sur)
summary(cox)
df <- with(sur, data.frame(pathologic_stage = c("Stage I", "Stage II", "Stage III", "Stage IV")))

cox <- coxph(Surv(dfs_days, dfs_status)~tils+pathologic_stage, data=sur)
summary(cox)
df <- with(sur, data.frame(tils = c("<20%", ">20%"), pathologic_stage = rep(levels(pathologic_stage)[2], 2)))

cox <- coxph(Surv(dfs_days, dfs_status)~tils+age_at_initial_pathologic_diagnosis+pathologic_stage, data=sur)
summary(cox)
#df <- with(sur, data.frame(tils = c("<1%", "1-20%", "20-40%", ">40%"), pathologic_stage = rep(levels(pathologic_stage)[2], 2)))
df <- with(sur, data.frame(tils = c("<20%", ">20%"), age_at_initial_pathologic_diagnosis = rep(mean(age_at_initial_pathologic_diagnosis), 2), 
                           pathologic_stage = rep(levels(pathologic_stage)[2], 2)))

summary(survfit(cox, newdata=df))
plot(survfit(cox, newdata=df))

width <- 600
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png("~/brca/analysis/R/dfs_coxph_tils.stage_sII.png", width=width, height=height)

palette <- c("#F8766D", "#00BFC4")
plot(survfit(cox, newdata=df), lty=c(1,1), main=paste("Disease Free Survival"),
     xlab="days", ylab="survival", col=palette, lwd=2, mark=108, cex=1.5)
legend(500,0.3, legend=c(paste0(">20% TILs"), paste0("<20% TILs")), 
       lty=c(1,1), col=rev(palette), lwd=2)

dev.off()


ggsurv(smpl.fit, xlab="days", size.est=0.75, cens.size=4, cens.shape=3)

#autoplot(smpl.fit)
