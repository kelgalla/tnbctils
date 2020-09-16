#---------------------------------------------------------------------
# FILE     : tnbctilsSurvMETABRIC.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-08
# COMMENTS : read in survival info and cibersort info and generate 
#            histograms and KM curves at different cut offs
#---------------------------------------------------------------------

library(ggplot2)
library(survival)
library(GGally)
library(ggfortify)

surFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcSurv.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
sur$OS_STATUS2 <- ifelse(sur$OS_STATUS == "0:LIVING", 0, 1)

vals <- "abs"
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

#-------------------------
# code for getting cibersort values for creating supp table 2
#round2 = function(x, n) {
#  posneg = sign(x)
#  z = abs(x)*10^n
#  z = z + 0.5 + sqrt(.Machine$double.eps)
#  z = trunc(z)
#  z = z/10^n
#  z*posneg
#}

#rel$T.cells.CD8.round <- round2(rel$T.cells.CD8, 3)
#rel$T.cells.CD4.memory.activated.round <- round2(rel$T.cells.CD4.memory.activated, 3)
#tmp <- rel[order(rel$Input.Sample), c(1,27,5,28,8)]
#----------------------

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\intensity\\"
dirSave2 <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\intensity\\noSig\\"

pvalTbl <- data.frame(cellType=character(), survival=character(), splitType=character(), lowEnd=numeric(), 
                      highEnd=numeric(), mean=numeric(), median=numeric(), sd=numeric(), pval=numeric())

# 5 = CD8, 8 = CD4
#cell_idx <- 11
for (cell_idx in 2:(dim(rel)[2]-3)){

  cibersort <- rel[,c(1,cell_idx)]
  cellType <- colnames(cibersort)[2]
  colnames(cibersort)[2] <- "sum"

  surv <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
  
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
  #times <- 0.25
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
    
    if (lowEnd <= 0){ 
      pvalTbl[nrow(pvalTbl)+1,] <- comp1
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

    #tmp <- surv[order(surv$PATIENT_ID), c(1,34,35)] 
    
    # km
    smpl.surv <- Surv(surv$OS_MONTHS, surv$OS_STATUS2)~surv$imm
    sigTest <- survdiff(smpl.surv)
    pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
    smpl.fit <- survfit(smpl.surv)
    
    comp1[9] <- pvalue
    pvalTbl[nrow(pvalTbl)+1,] <- comp1
    
    ifelse(pvalue < 0.05, save <- dirSave, save <- dirSave2)
      
    #---
    # os
    width <- 600
    #ratio <- 461/572
    ratio <- 547/600
    height <- trunc(width*ratio)
    png(paste0(save, "os.", vals, ".", cellType, ".", times, 
               ifelse(times=="mean" | times=="median", "", "sd"), ".png"), width=width, height=height)
      
    palette <- c("#F8766D", "#00BFC4")
    plot(smpl.fit, lty=c(1,1), main=paste0("Overall Survival ", cellType, "\nPvalue =", round(pvalue, digits=3)),
         xlab="months", ylab="survival", col=palette, lwd=2, mark=108, cex=1.5)
    legend(50,0.3, legend=c(paste0(">", toString(highEnd), " ", cellType, " n=", highNum), paste0("<", toString(lowEnd), " ", cellType, " n=", lowNum)), 
           lty=c(1,1), col=rev(palette), lwd=2)
      
    dev.off()
      
    #---
    # fancy plot
    # os
    width <- 2500
    #ratio <- 461/572
    ratio <- 547/600
    height <- trunc(width*ratio)
    png(paste0(save, "os.", vals, ".", cellType, ".", times, 
               ifelse(times=="mean" | times=="median", "", "sd"), ".2.png"), 
        width=width, height=height, res=300)
      
    palette <- c("#F8766D", "#00BFC4")
    par(mar=c(5.1,5.1,4.1,2.1))
    plot(smpl.fit, lty=c(1,1),
         xlab="Months", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
    legend(50,0.35, legend=c(paste0("High, n=", highNum), paste0("Low, n=", lowNum)), 
           lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)
      
    dev.off()
    # fancy plot
    #---
    
  }
}

cellTblOS <- pvalTbl[pvalTbl$survival=="os" & pvalTbl$splitType=="0.25",]
cellTblOS$fdr <- p.adjust(cellTblOS$pval, method="fdr")


