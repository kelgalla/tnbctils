#---------------------------------------------------------------------
# FILE     : clusterGrp.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2019-07-21
# COMMENTS : read in a tab-delimited file with TCGA info, including 
#            survival and cluster information and generate KM survival 
#            curves
#---------------------------------------------------------------------

library(survival)

surFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))

clusterGrpFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\clusterGrp.txt"
clusterGrp <- read.delim(clusterGrpFile, header=T, stringsAsFactors=F)
clusterGrp$barcode <- gsub("\\.", "-", clusterGrp$barcode)

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\fpkm\\"

sur <- merge(sur, clusterGrp, by.x="bcr_patient_barcode", by.y="barcode", all.x=T)
sur <- sur[!(is.na(sur$clusterGrp)),]

sur <- sur[sur$clusterGrp %in% c("bothNeg", "bothPos"),]
#sur <- sur[sur$clusterGrp %in% c("bothNeg", "cd8upcd4down"),]
#sur <- sur[sur$clusterGrp %in% c("bothNeg", "cd4upcd8down"),]
#sur$clusterGrp <- factor(sur$clusterGrp, levels=c("bothNeg", "bothPos", "cd8upcd4down", "cd4upcd8down"))

# bothPos #EB00FF
# cd8upcd4down #00FF0A
# cd4upcd8down #FF0099
# bothNeg #FF9900

smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$clusterGrp
sigTest <- survdiff(smpl.surv)
sigTest
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)

highNum <- length(sur[sur$clusterGrp=="bothPos",]$clusterGrp)
lowNum <- length(sur[sur$clusterGrp=="bothNeg",]$clusterGrp)

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "os.bothPosvsbothNeg.png"), 
    width=width, height=height, res=300)

palette <- c("#FF9900FF", "#EB00FFFF")
par(mar=c(5.1,5.1,4.1,2.1))
plot(smpl.fit, lty=c(1,1),
     xlab="Days", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(450,0.35, legend=c(paste0("CD8/CD4 High, n=", highNum), paste0("CD8/CD4 Low, n=", lowNum)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()

cox <- coxph(Surv(os_days, os_status)~clusterGrp, data=sur)
summary(cox)
