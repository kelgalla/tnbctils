#---------------------------------------------------------------------
# FILE     : tnbctilsInfoMETABRIC.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-15
# COMMENTS : read in a tab-delimited file with TCGA info, including 
#            survival information, and summarize it for table 1, also 
#            univariate analysis for cell specific
#---------------------------------------------------------------------

library(survival)
library(plyr)
library(survminer)

surFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcSurv.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
sur$OS_STATUS2 <- ifelse(sur$OS_STATUS == "0:LIVING", 0, 1)

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\intensity\\"
dirSave2 <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\intensity\\noSig\\"

#----------------------------------
# 199 samples
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
#sur <- sur[!(is.na(sur$sum)),] # want data on 199 and not 196

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

mean(sur$AGE_AT_DIAGNOSIS)

#----------------------------------
# CD8 continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

# cox os
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum, data=sur)
summary(cox)

cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum+TUMOR_STAGE+AGE_AT_DIAGNOSIS+INFERRED_MENOPAUSAL_STATE, data=sur)
summary(cox)

#----------------------------------
# CD8 quartiles
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

quantile(sur$sum)

sur$quartile <- cut(sur$sum, breaks=quantile(sur$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)

#tmp <- sur[order(sur$sum),c(1,34,35)]

q1Num <- dim(sur[sur$quartile=="Q1",])[1]
q2Num <- dim(sur[sur$quartile=="Q2",])[1]
q3Num <- dim(sur[sur$quartile=="Q3",])[1]
q4Num <- dim(sur[sur$quartile=="Q4",])[1]

# km OS
smpl.surv <- Surv(sur$OS_MONTHS, sur$OS_STATUS2)~sur$quartile
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

#---
# fancy plot
# os
width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave2, "os.abs.T.cells.CD8.quartiles.2.png"), 
    width=width, height=height, res=300)

palette <- c("#B79F00", "#00BA38", "#619CFF", "#A58AFF")
par(mar=c(5.1,5.1,4.1,2.1))
plot(smpl.fit, lty=c(1,1),
     xlab="Months", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(200, 1.03, legend=c(paste0("Q4, n=", q1Num), paste0("Q3, n=", q2Num), paste0("Q2, n=", q3Num), paste0("Q1, n=", q4Num)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()
# fancy plot
#---

#n <- 7
#hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)

# cox os
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~quartile, data=sur)
summary(cox)

#----------------------------------
# CD8 high and low
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]
mean(sur$AGE_AT_DIAGNOSIS) # for 196

m <- mean(sur$sum, na.rm=T)
sd <- sd(sur$sum, na.rm=T)
times <- 0.25
lowEnd <- round(m - (times * sd), digits=4)
highEnd <- round(m + (times * sd), digits=4)

# lowEnd, highEnd
sur$imm <- NA
CDgreaterThan <- paste0(">", toString(highEnd))
sur[sur$sum > highEnd & !is.na(sur$sum),]$imm <- CDgreaterThan
CDlessThan <- paste0("<", toString(lowEnd))
sur[sur$sum < lowEnd & !is.na(sur$sum),]$imm <- CDlessThan
sur <- sur[!(is.na(sur$imm)),]
sur$imm <- factor(sur$imm, levels=c(paste0("<", toString(lowEnd)), paste0(">", toString(highEnd))), ordered=T)
lowNum <- dim(sur[sur$imm==paste0("<", toString(lowEnd)),])[1]
highNum <- dim(sur[sur$imm==paste0(">", toString(highEnd)),])[1]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

greaterThan <- sur[sur$imm == CDgreaterThan,]
mean(greaterThan$AGE_AT_DIAGNOSIS)
lessThan <- sur[sur$imm == CDlessThan,]
mean(lessThan$AGE_AT_DIAGNOSIS)
wilcox.test(greaterThan$AGE_AT_DIAGNOSIS, 
            lessThan$AGE_AT_DIAGNOSIS)

count(greaterThan, 'TUMOR_STAGE')
count(lessThan, 'TUMOR_STAGE')
t <- matrix(c(42,4,19,67,7,27), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("high CD8", "low CD8")))
fisher.test(t)

count(greaterThan, 'INFERRED_MENOPAUSAL_STATE')
count(lessThan, 'INFERRED_MENOPAUSAL_STATE')
t <- matrix(c(28,37,42,59), 2, 2, dimnames=list(c("pre", "post"), c("high CD8", "low CD8")))
fisher.test(t)

# km OS
smpl.surv <- Surv(sur$OS_MONTHS, sur$OS_STATUS2)~sur$imm
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
pvalue

summary(smpl.fit, times=60) # 5 year survival
summary(smpl.fit, times=84) # 7 year survival
summary(smpl.fit, times=120) # 10 year survival
summary(smpl.fit, times=180) # 15 year survival
summary(smpl.fit, times=240) # 20 year survival

# cox OS
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~imm, data=sur)
summary(cox)

#----------------------------------
# T Cells gamma delta continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 11
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

# cox os
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum, data=sur)
summary(cox)

cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum+TUMOR_STAGE+AGE_AT_DIAGNOSIS+INFERRED_MENOPAUSAL_STATE, data=sur)
summary(cox)

cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum+AGE_AT_DIAGNOSIS, data=sur)
summary(cox)


#----------------------------------
# T Cells gamma delta quartiles
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 11
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

quantile(sur$sum)

sur$quartile <- cut(sur$sum, breaks=quantile(sur$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)

#tmp <- sur[order(sur$sum),c(1,34,35)]

q1Num <- dim(sur[sur$quartile=="Q1",])[1]
q2Num <- dim(sur[sur$quartile=="Q2",])[1]
q3Num <- dim(sur[sur$quartile=="Q3",])[1]
q4Num <- dim(sur[sur$quartile=="Q4",])[1]

# km OS
smpl.surv <- Surv(sur$OS_MONTHS, sur$OS_STATUS2)~sur$quartile
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

#---
# fancy plot
# os
width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "os.abs.T.cells.gamma.delta.quartiles.2.png"), 
    width=width, height=height, res=300)

palette <- c("#B79F00", "#00BA38", "#619CFF", "#A58AFF")
par(mar=c(5.1,5.1,4.1,2.1))
par(xpd=TRUE)
plot(smpl.fit, lty=c(1,1),
     xlab="Months", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(-10, 0.38, legend=c(paste0("Q4, n=", q1Num), paste0("Q3, n=", q2Num), paste0("Q2, n=", q3Num), paste0("Q1, n=", q4Num)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()
# fancy plot
#---

#n <- 7
#hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)

# cox os
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~quartile, data=sur)
summary(cox)

#----------------------------------
# T Cells gamma delta high and low
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 11
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]
mean(sur$AGE_AT_DIAGNOSIS) # for 196

m <- mean(sur$sum, na.rm=T)
sd <- sd(sur$sum, na.rm=T)
times <- 0.25
lowEnd <- round(m - (times * sd), digits=4)
highEnd <- round(m + (times * sd), digits=4)

# lowEnd, highEnd
sur$imm <- NA
CDgreaterThan <- paste0(">", toString(highEnd))
sur[sur$sum > highEnd & !is.na(sur$sum),]$imm <- CDgreaterThan
CDlessThan <- paste0("<", toString(lowEnd))
sur[sur$sum < lowEnd & !is.na(sur$sum),]$imm <- CDlessThan
sur <- sur[!(is.na(sur$imm)),]
sur$imm <- factor(sur$imm, levels=c(paste0("<", toString(lowEnd)), paste0(">", toString(highEnd))), ordered=T)
lowNum <- dim(sur[sur$imm==paste0("<", toString(lowEnd)),])[1]
highNum <- dim(sur[sur$imm==paste0(">", toString(highEnd)),])[1]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

greaterThan <- sur[sur$imm == CDgreaterThan,]
mean(greaterThan$AGE_AT_DIAGNOSIS)
lessThan <- sur[sur$imm == CDlessThan,]
mean(lessThan$AGE_AT_DIAGNOSIS)
wilcox.test(greaterThan$AGE_AT_DIAGNOSIS, 
            lessThan$AGE_AT_DIAGNOSIS)

count(greaterThan, 'TUMOR_STAGE')
count(lessThan, 'TUMOR_STAGE')
t <- matrix(c(44,4,18,65,6,32), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("high CD8", "low CD8")))
fisher.test(t)

count(greaterThan, 'INFERRED_MENOPAUSAL_STATE')
count(lessThan, 'INFERRED_MENOPAUSAL_STATE')
t <- matrix(c(27,39,43,60), 2, 2, dimnames=list(c("pre", "post"), c("high CD8", "low CD8")))
fisher.test(t)

# km OS
smpl.surv <- Surv(sur$OS_MONTHS, sur$OS_STATUS2)~sur$imm
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
pvalue

summary(smpl.fit, times=60) # 5 year survival
summary(smpl.fit, times=84) # 7 year survival
summary(smpl.fit, times=120) # 10 year survival
summary(smpl.fit, times=180) # 15 year survival
summary(smpl.fit, times=240) # 20 year survival

# cox OS
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~imm, data=sur)
summary(cox)

cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~imm+AGE_AT_DIAGNOSIS, data=sur)
summary(cox)


#----------------------------------
# CD4 continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

# cox os
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum, data=sur)
summary(cox)

cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum+TUMOR_STAGE+AGE_AT_DIAGNOSIS+INFERRED_MENOPAUSAL_STATE, data=sur)
summary(cox)

#----------------------------------
# CD4 high and low
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]
mean(sur$AGE_AT_DIAGNOSIS) # for 196

m <- mean(sur$sum, na.rm=T)
sd <- sd(sur$sum, na.rm=T)
times <- 0.25
lowEnd <- round(m - (times * sd), digits=4)
highEnd <- round(m + (times * sd), digits=4)

# lowEnd, highEnd
sur$imm <- NA
CDgreaterThan <- paste0(">", toString(highEnd))
sur[sur$sum > highEnd & !is.na(sur$sum),]$imm <- CDgreaterThan
CDlessThan <- paste0("<", toString(lowEnd))
sur[sur$sum < lowEnd & !is.na(sur$sum),]$imm <- CDlessThan
sur <- sur[!(is.na(sur$imm)),]
sur$imm <- factor(sur$imm, levels=c(paste0("<", toString(lowEnd)), paste0(">", toString(highEnd))), ordered=T)
lowNum <- dim(sur[sur$imm==paste0("<", toString(lowEnd)),])[1]
highNum <- dim(sur[sur$imm==paste0(">", toString(highEnd)),])[1]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

greaterThan <- sur[sur$imm == CDgreaterThan,]
mean(greaterThan$AGE_AT_DIAGNOSIS)
lessThan <- sur[sur$imm == CDlessThan,]
mean(lessThan$AGE_AT_DIAGNOSIS)
wilcox.test(greaterThan$AGE_AT_DIAGNOSIS, 
            lessThan$AGE_AT_DIAGNOSIS)

count(greaterThan, 'TUMOR_STAGE')
count(lessThan, 'TUMOR_STAGE')
t <- matrix(c(37,0,16,77,10,34), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("high CD4", "low CD4")))
fisher.test(t)

count(greaterThan, 'INFERRED_MENOPAUSAL_STATE')
count(lessThan, 'INFERRED_MENOPAUSAL_STATE')
t <- matrix(c(23,30,47,74), 2, 2, dimnames=list(c("pre", "post"), c("high CD4", "low CD4")))
fisher.test(t)

# cox OS
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~imm, data=sur)
summary(cox)

#----------------------------------
# M1 macrophage continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 16
cibersort <- rel[,c(1,cell_idx)]
cellType <- colnames(cibersort)[2]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$TUMOR_STAGE %in% c("1", "2"),]$TUMOR_STAGE <- "Stage I-II"
sur[sur$TUMOR_STAGE %in% c("3"),]$TUMOR_STAGE <- "Stage III-IV"
sur$TUMOR_STAGE <- factor(sur$TUMOR_STAGE, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'TUMOR_STAGE')

sur$INFERRED_MENOPAUSAL_STATE <- factor(sur$INFERRED_MENOPAUSAL_STATE, levels=c("Pre", "Post"), ordered=T)
count(sur, 'INFERRED_MENOPAUSAL_STATE')

# cox os
cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum, data=sur)
summary(cox)

#----------------------------------
# all immune types continuous for purposes of calculating fdr
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

pvalTbl <- data.frame(cellType=character(), survival=character(), splitType=character(), pval=numeric())

for (cell_idx in 2:(dim(rel)[2]-3)){
  
  cibersort <- rel[,c(1,cell_idx)]
  cellType <- colnames(cibersort)[2]
  colnames(cibersort)[2] <- "sum"
  
  surv <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
  
  # analysis of cibersort
  surv <- surv[!(is.na(surv$sum)),]
  
  # cox os
  cox <- coxph(Surv(OS_MONTHS, OS_STATUS2)~sum, data=surv)
  summary(cox)
  pvalue <- summary(cox)$coefficients[,5]
  
  comp1 <- c(cellType, "os", "continuous", pvalue)
  pvalTbl[nrow(pvalTbl)+1,] <- comp1
  
}

cellTblOS <- pvalTbl[pvalTbl$survival=="os",]
cellTblOS$fdr <- p.adjust(cellTblOS$pval, method="fdr")

#----------------------------------
# all immune types quartiles for purposes of calculating fdr
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

pvalTbl <- data.frame(cellType=character(), survival=character(), splitType=character(), pval=numeric())

#cell_idx=2
for (cell_idx in 2:(dim(rel)[2]-3)){
  
  cibersort <- rel[,c(1,cell_idx)]
  cellType <- colnames(cibersort)[2]
  colnames(cibersort)[2] <- "sum"
  
  surv <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
  
  # analysis of cibersort
  surv <- surv[!(is.na(surv$sum)),]
  
  p <- tryCatch({ # some immune cells have too many 0's such that they can't be split into quartiles
    
    surv$quartile <- cut(surv$sum, breaks=quantile(surv$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)
    
    # km OS
    smpl.surv <- Surv(surv$OS_MONTHS, surv$OS_STATUS2)~surv$quartile
    sigTest <- survdiff(smpl.surv)
    pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
    
  }, error = function(e) {
    
    return(NA)
    
  })
  
  comp1 <- c(cellType, "os", "quartile", p)
  pvalTbl[nrow(pvalTbl)+1,] <- comp1
  
}

cellTblOS <- pvalTbl[pvalTbl$survival=="os",]
cellTblOS$fdr <- p.adjust(cellTblOS$pval, method="fdr")



