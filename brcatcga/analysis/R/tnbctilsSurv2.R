#---------------------------------------------------------------------
# FILE     : tnbctilsSurv2.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-08-05
# COMMENTS : read in survival info and calculate univariate and 
#            generate KM curves for tils analysis
#---------------------------------------------------------------------

library(survival)

surFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
sur <- sur[!(is.na(sur$tils)),]

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\"

#---------------------
# analysis of tils
#----------------------
#sur$tils <- factor(sur$tils, levels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"), ordered=T)

# continuous
sur[sur$tils %in% c("<1%"),]$tils <- 0
sur[sur$tils %in% c("1-10%"),]$tils <- 1
sur[sur$tils %in% c("10-20%"),]$tils <- 2
sur[sur$tils %in% c("20-30%"),]$tils <- 3
sur[sur$tils %in% c("30-40%"),]$tils <- 4
sur[sur$tils %in% c("40-50%"),]$tils <- 5
sur[sur$tils %in% c("50-60%"),]$tils <- 6
sur[sur$tils %in% c(">70%"),]$tils <- 7
sur$tils <- as.numeric(sur$tils)

# cox
cox <- coxph(Surv(os_days, os_status)~tils, data=sur)
summary(cox)

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~tils, data=sur)
summary(cox)

# <30% and >30%
sur[sur$tils %in% c("<1%", "1-10%", "10-20%", "20-30%"),]$tils <- "<30%"
sur[sur$tils %in% c("30-40%", "40-50%", "50-60%", "60-70%", ">70%"),]$tils <- ">30%"
sur$tils <- factor(sur$tils, levels=c("<30%", ">30%"), ordered=T)
lowNum <- dim(sur[sur$tils=="<30%",])[1]
highNum <- dim(sur[sur$tils==">30%",])[1]

# km os
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$tils
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "os2-30.png"), width=width, height=height, res=300)


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
png(paste0(dirSave, "dfs2-30.png"), 
    width=width, height=height, res=300)

palette <- c("#F8766D", "#00BFC4")
par(mar=c(5.1,5.1,4.1,2.1))
plot(smpl.fit, lty=c(1,1),
     xlab="Days", ylab="Disease Free Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(500,0.35, legend=c(paste0(">30% TILs, n=", highNum), paste0("<30% TILs, n=", lowNum)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()






