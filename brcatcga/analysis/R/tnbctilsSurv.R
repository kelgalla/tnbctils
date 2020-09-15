#---------------------------------------------------------------------
# FILE     : tnbctilsSurv.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-08-05
# COMMENTS : read in survival info and cibersort info and generate 
#            histograms and KM curves at different cut offs
#---------------------------------------------------------------------

library(ggplot2)
library(survival)
library(GGally)
library(ggfortify)

surFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))

vals <- "abs"
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\fpkm\\"

pvalTbl <- data.frame(cellType=character(), survival=character(), splitType=character(), lowEnd=numeric(), 
                      highEnd=numeric(), mean=numeric(), median=numeric(), sd=numeric(), pval=numeric())

# 5 = CD8, 8 = CD4
#cell_idx <- 2
for (cell_idx in 2:(dim(rel)[2]-3)){
  
  cibersort <- rel[,c(1,cell_idx)]
  cellType <- colnames(cibersort)[2]
  colnames(cibersort)[2] <- "sum"
  
  surv <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
  
  # analysis of cibersort
  surv <- surv[!(is.na(surv$sum)),]
  
  x <- surv$sum
  
  width <- 600
  #ratio <- 461/572
  ratio <- 547/600
  height <- trunc(width*ratio)
  png(paste0(dirSave, "hist.", vals, ".", cellType, ".png"), width=width, height=height)
  
  # automated histo
  h<-hist(x, breaks=10, col="red", xlab="quantity of cell type",
          main=paste0("Histogram of ", cellType))
  xfit<-seq(min(x),max(x),length=40)
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
  yfit <- yfit*diff(h$mids[1:2])*length(x)
  lines(xfit, yfit, col="blue", lwd=2)
  
  dev.off()
  
  # mean +- times * sd
  m <- mean(surv$sum, na.rm=T)
  med <- median(surv$sum, na.rm=T)
  sd <- sd(surv$sum, na.rm=T)
#  times <- 0.25
  for (times in c(0.25, 0.5, 1, "mean", "median")){
    #print(paste0(cellType, " ", times))
    lowEnd <- -1
    highEnd <- -1
    if (times == "mean"){
      lowEnd <- m
      highEnd <- m
    } else if (times == "median"){
      lowEnd <- med
      highEnd <- med
    } else {
      lowEnd <- round(m - (as.numeric(times) * sd), digits=4)
      highEnd <- round(m + (as.numeric(times) * sd), digits=4)
    }
    
    comp1 <- c(cellType, "os", times, lowEnd, highEnd, m, med, sd, NA)
    comp2 <- c(cellType, "dfs", times, lowEnd, highEnd, m, med, sd, NA)
    
    if (lowEnd <= 0){ 
      pvalTbl[nrow(pvalTbl)+1,] <- comp1
      pvalTbl[nrow(pvalTbl)+1,] <- comp2
      next
      }
    
    # lowEnd, highEnd
    surv$imm <- NA
    surv[surv$sum < lowEnd & !is.na(surv$sum),]$imm <- paste0("<", toString(lowEnd))
    surv[surv$sum > highEnd & !is.na(surv$sum),]$imm <- paste0(">", toString(highEnd))
    surv <- surv[!(is.na(surv$imm)),]
    surv$imm <- factor(surv$imm, levels=c(paste0("<", toString(lowEnd)), paste0(">", toString(highEnd))), ordered=T)
    lowNum <- dim(surv[surv$imm==paste0("<", toString(lowEnd)),])[1]
    highNum <- dim(surv[surv$imm==paste0(">", toString(highEnd)),])[1]
    
    # km
    smpl.surv <- Surv(surv$os_days, surv$os_status)~surv$imm
    sigTest <- survdiff(smpl.surv)
    pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
    smpl.fit <- survfit(smpl.surv)
    
    comp1[9] <- pvalue
    pvalTbl[nrow(pvalTbl)+1,] <- comp1

    #if (pvalue < 0.05) {
      
      #---
      # os
      width <- 600
      #ratio <- 461/572
      ratio <- 547/600
      height <- trunc(width*ratio)
      png(paste0(dirSave, "os.", vals, ".", cellType, ".", times, 
                 ifelse(times=="mean" | times=="median", "", "sd"), ".png"), width=width, height=height)
      
      palette <- c("#F8766D", "#00BFC4")
      plot(smpl.fit, lty=c(1,1), main=paste0("Overall Survival ", cellType, "\nPvalue =", round(pvalue, digits=3)),
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
      png(paste0(dirSave, "os.", vals, ".", cellType, ".", times, 
                 ifelse(times=="mean" | times=="median", "", "sd"), ".2.png"), 
          width=width, height=height, res=300)
      
      palette <- c("#F8766D", "#00BFC4")
      par(mar=c(5.1,5.1,4.1,2.1))
      plot(smpl.fit, lty=c(1,1),
           xlab="Days", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
      legend(500,0.35, legend=c(paste0("High, n=", highNum), paste0("Low, n=", lowNum)), 
             lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)
      
      dev.off()
      # fancy plot
      #---
    #}
    
    # km
    smpl.surv <- Surv(surv$dfs_days, surv$dfs_status)~surv$imm
    sigTest <- survdiff(smpl.surv)
    pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
    smpl.fit <- survfit(smpl.surv)
    
    comp2[9] <- pvalue
    pvalTbl[nrow(pvalTbl)+1,] <- comp2
    
    #if (pvalue < 0.05){
      
      #---
      # dfs
      width <- 600
      #ratio <- 461/572
      ratio <- 547/600
      height <- trunc(width*ratio)
      png(paste0(dirSave, "dfs.", vals, ".", cellType, ".", times, 
                 ifelse(times=="mean" | times=="median", "", "sd"), ".png"), width=width, height=height)
      
      palette <- c("#F8766D", "#00BFC4")
      plot(smpl.fit, lty=c(1,1), main=paste0("Disease Free Survival ", cellType, "\nPvalue =", round(pvalue, digits=3)),
           xlab="days", ylab="survival", col=palette, lwd=2, mark=108, cex=1.5)
      legend(500,0.3, legend=c(paste0(">", toString(highEnd), " ", cellType, " n=", highNum), paste0("<", toString(lowEnd), " ", cellType, " n=", lowNum)), 
             lty=c(1,1), col=rev(palette), lwd=2)
      
      dev.off()
      
      #---
      # fancy plot
      # dfs
      width <- 2500
      #ratio <- 461/572
      ratio <- 547/600
      height <- trunc(width*ratio)
      png(paste0(dirSave, "dfs.", vals, ".", cellType, ".", times,
                 ifelse(times=="mean" | times=="median", "", "sd"), ".2.png"), 
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
      
    }
  }

cellTblOS <- pvalTbl[pvalTbl$survival=="os" & pvalTbl$splitType=="0.25",]
cellTblOS$fdr <- p.adjust(cellTblOS$pval, method="fdr")

cellTblDFS <- pvalTbl[pvalTbl$survival=="dfs" & pvalTbl$splitType=="0.25",]
cellTblDFS$fdr <- p.adjust(cellTblDFS$pval, method="fdr")
