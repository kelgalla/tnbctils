#---------------------------------------------------------------------
# FILE     : tnbctilsInfo.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2019-12-20
# COMMENTS : read in a tab-delimited file with TCGA info, including 
#            survival information, and summarize it for table 1, also 
#            univariate and multivariate analysis for cell specific
#---------------------------------------------------------------------

library(survival)
library(plyr)
library(dplyr)
library(survminer)

surFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\fpkm\\"

#----------------------------------
# 133 samples
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
#sur <- sur[!(is.na(sur$sum)),] # want data on 133 and not 132

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

mean(sur$age_at_initial_pathologic_diagnosis)

#----------------------------------
# tils > and < 30%
#----------------------------------
sur <- sur[!(is.na(sur$tils)),]

sur[sur$tils %in% c("<1%", "1-10%", "10-20%", "20-30%"),]$tils <- "<30%"
sur[sur$tils %in% c("30-40%", "40-50%", "50-60%", "60-70%", ">70%"),]$tils <- ">30%"
sur$tils <- factor(sur$tils, levels=c("<30%", ">30%"), ordered=T)
lowNum <- dim(sur[sur$tils=="<30%",])[1]
highNum <- dim(sur[sur$tils==">30%",])[1]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

greaterThan30 <- sur[sur$tils==">30%",]
mean(greaterThan30$age_at_initial_pathologic_diagnosis)
lessThan30 <- sur[sur$tils=="<30%",]
mean(lessThan30$age_at_initial_pathologic_diagnosis)
wilcox.test(greaterThan30$age_at_initial_pathologic_diagnosis, lessThan30$age_at_initial_pathologic_diagnosis)
#cor.test(greaterThan30$age_at_initial_pathologic_diagnosis, lessThan30$age_at_initial_pathologic_diagnosis)

count(greaterThan30, 'race')
count(lessThan30, 'race')

t <- matrix(c(6,4,1,52,30,10), 3, 2, dimnames=list(c("white", "black", "na"), c("tils>30", "tils<30")))
#t <- matrix(c(6,4,1,52,30,10), 2, 3, byrow=T, dimnames=list(c("tils>30", "tils<30"), c("white", "black", "na")))
fisher.test(t)

count(greaterThan30, 'pathologic_stage')
count(lessThan30, 'pathologic_stage')
t <- matrix(c(11, 0, 0, 73, 16, 3), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("tils>30", "tils<30")))
fisher.test(t)

count(greaterThan30, 'menopause_status')
count(lessThan30, 'menopause_status')
t <- matrix(c(4,5,2,21,59,12), 3, 2, dimnames=list(c("pre", "post", "na"), c("tils>30", "tils<30")))
fisher.test(t)

# km OS
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$tils
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
pvalue

# cox OS
cox <- coxph(Surv(os_days, os_status)~tils, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

# km DFS
smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$tils
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
pvalue

# cox DFS
cox <- coxph(Surv(dfs_days, dfs_status)~tils, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

#----------------------------------
# CD8 continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

# cox os
cox <- coxph(Surv(os_days, os_status)~sum, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~sum+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~sum+pathologic_stage+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~sum, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~sum+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~sum+pathologic_stage, data=sur)
summary(cox)

#----------------------------------
# CD8 quartiles
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

quantile(sur$sum)

sur$quartile <- cut(sur$sum, breaks=quantile(sur$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)

#tmp <- sur[order(sur$sum),c(1,75,76)]

q1Num <- dim(sur[sur$quartile=="Q1",])[1]
q2Num <- dim(sur[sur$quartile=="Q2",])[1]
q3Num <- dim(sur[sur$quartile=="Q3",])[1]
q4Num <- dim(sur[sur$quartile=="Q4",])[1]

#sur$quartileNum <- ntile(sur$sum, 4)
#sur$quartile <- NA
#sur[sur$quartileNum==1,]$quartile <- "Q1"
#sur[sur$quartileNum==2,]$quartile <- "Q2"
#sur[sur$quartileNum==3,]$quartile <- "Q3"
#sur[sur$quartileNum==4,]$quartile <- "Q4"
#sur$quartile <- factor(sur$quartile, levels=c("Q1", "Q2", "Q3", "Q4"))

# km OS
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$quartile
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
png(paste0(dirSave, "os.abs.T.cells.CD8.quartiles.2.png"), 
    width=width, height=height, res=300)

palette <- c("#B79F00", "#00BA38", "#619CFF", "#A58AFF")
par(mar=c(5.1,5.1,4.1,2.1))
plot(smpl.fit, lty=c(1,1),
     xlab="Days", ylab="Overall Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(500,0.5, legend=c(paste0("Q4, n=", q1Num), paste0("Q3, n=", q2Num), paste0("Q2, n=", q3Num), paste0("Q1, n=", q4Num)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()
# fancy plot
#---

#n <- 7
#hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)

# cox os
cox <- coxph(Surv(os_days, os_status)~quartile, data=sur)
summary(cox)

# km DFS
smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$quartile
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
survfit(smpl.surv)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~quartile, data=sur)
summary(cox)


#----------------------------------
# CD8 high and low
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]
mean(sur$age_at_initial_pathologic_diagnosis) # for 132

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

#sur$raceEth <- with(sur, paste0(sur$race, " | ", sur$ethnicity))
#sur$raceEth <- factor(sur$raceEth, levels=c(unique(sur$raceEth)))
#unique(sur$raceEth)

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

#sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- "WHITE OR ASIAN"
#sur[sur$race=="WHITE" & !is.na(sur$race),]$race <- "WHITE OR ASIAN"
#sur$race <- factor(sur$race, levels=c(unique(sur$race)))
#count(sur, 'race')

#sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I"
#sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage II"
#sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III"
#sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I", "Stage II", "Stage III"), ordered=T)
#count(sur, 'pathologic_stage')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

greaterThan <- sur[sur$imm == CDgreaterThan,]
mean(greaterThan$age_at_initial_pathologic_diagnosis)
lessThan <- sur[sur$imm == CDlessThan,]
mean(lessThan$age_at_initial_pathologic_diagnosis)
wilcox.test(greaterThan$age_at_initial_pathologic_diagnosis, 
            lessThan$age_at_initial_pathologic_diagnosis)

count(greaterThan, 'race')
count(lessThan, 'race')
t <- matrix(c(24,14,2,37,27,7), 3, 2, dimnames=list(c("white", "black", "na"), c("high CD8", "low CD8")))
fisher.test(t)

count(greaterThan, 'pathologic_stage')
count(lessThan, 'pathologic_stage')
t <- matrix(c(35,5,0,57,12,2), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("high CD8", "low CD8")))
fisher.test(t)

count(greaterThan, 'menopause_status')
count(lessThan, 'menopause_status')
t <- matrix(c(12,24,4,18,43,10), 3, 2, dimnames=list(c("pre", "post", "na"), c("high CD8", "low CD8")))
fisher.test(t)

# km OS
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$imm
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
survfit(smpl.surv)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

summary(smpl.fit, times=1826.25) # 5 year survival
summary(smpl.fit, times=2556.75) # 7 year survival
summary(smpl.fit, times=3652.5, extend=T) # 10 year survival

# cox OS
cox <- coxph(Surv(os_days, os_status)~imm, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~race, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~pathologic_stage, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~menopause_status, data=sur)
summary(cox)

# all
cox <- coxph(Surv(os_days, os_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~imm+pathologic_stage+race+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

prop <- cox.zph(cox, transform="km", global=T)

print(prop)
plot(prop)

ggforest(cox, data=sur)

# km DFS
#smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$imm
#sigTest <- survdiff(smpl.surv)
#pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
#smpl.fit <- survfit(smpl.surv)
#pvalue

# cox DFS
cox <- coxph(Surv(dfs_days, dfs_status)~imm, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

#----------------------------------
# CD4 continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

# cox os
cox <- coxph(Surv(os_days, os_status)~sum, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~sum+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~sum, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~sum+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

#----------------------------------
# CD4 quartiles
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

quantile(sur$sum)

sur$quartile <- cut(sur$sum, breaks=quantile(sur$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)

q1Num <- dim(sur[sur$quartile=="Q1",])[1]
q2Num <- dim(sur[sur$quartile=="Q2",])[1]
q3Num <- dim(sur[sur$quartile=="Q3",])[1]
q4Num <- dim(sur[sur$quartile=="Q4",])[1]

# km OS
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$quartile
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

#n <- 7
#hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)

# cox os
cox <- coxph(Surv(os_days, os_status)~quartile, data=sur)
summary(cox)

# km DFS
smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$quartile
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
survfit(smpl.surv)
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
png(paste0(dirSave, "dfs.abs.T.cells.CD4.memory.activated.quartiles.2.png"), 
    width=width, height=height, res=300)

palette <- c("#B79F00", "#00BA38", "#619CFF", "#A58AFF")
par(mar=c(5.1,5.1,4.1,2.1))
plot(smpl.fit, lty=c(1,1),
     xlab="Days", ylab="Disease Free Survival", col=palette, lwd=2, mark=108, cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
legend(200,0.5, legend=c(paste0("Q4, n=", q1Num), paste0("Q3, n=", q2Num), paste0("Q2, n=", q3Num), paste0("Q1, n=", q4Num)), 
       lty=c(1,1), col=rev(palette), lwd=2, cex=1.5)

dev.off()
# fancy plot
#---

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~quartile, data=sur)
summary(cox)

#----------------------------------
# CD4 high and low
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]
mean(sur$age_at_initial_pathologic_diagnosis) # for 132

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

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

#sur[sur$age_at_initial_pathologic_diagnosis > 30 & sur$age_at_initial_pathologic_diagnosis <= 40,]$age_at_initial_pathologic_diagnosis <- "<=40"
#sur[sur$age_at_initial_pathologic_diagnosis > 40,]$age_at_initial_pathologic_diagnosis <- ">40"
#sur[sur$age_at_initial_pathologic_diagnosis >= 70,]$age_at_initial_pathologic_diagnosis <- ">=70"
#sur$age_at_initial_pathologic_diagnosis <- factor(sur$age_at_initial_pathologic_diagnosis, levels=c("<=40", ">40"), ordered=T)
#count(sur, 'age_at_initial_pathologic_diagnosis')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

greaterThan <- sur[sur$imm == CDgreaterThan,]
mean(greaterThan$age_at_initial_pathologic_diagnosis)
lessThan <- sur[sur$imm == CDlessThan,]
mean(lessThan$age_at_initial_pathologic_diagnosis)
wilcox.test(greaterThan$age_at_initial_pathologic_diagnosis, 
            lessThan$age_at_initial_pathologic_diagnosis)

count(greaterThan, 'race')
count(lessThan, 'race')
t <- matrix(c(27,12,6,39,29,4), 3, 2, dimnames=list(c("white", "black", "na"), c("high CD4", "low CD4")))
fisher.test(t)

count(greaterThan, 'pathologic_stage')
count(lessThan, 'pathologic_stage')
pvals <- c()
t <- matrix(c(42,2,1,55,16,1), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("high CD4", "low CD4")))
fisher.test(t)
t <- matrix(c(42,2,55,16), 2, 2, dimnames=list(c("stage I-II", "stage III-IV"), c("high CD4", "low CD4")))
ftst <- fisher.test(t)
pvals <- c(pvals, ftst$p.value)
t <- matrix(c(2,1,16,1), 2, 2, dimnames=list(c("stage III-IV", "na"), c("high CD4", "low CD4")))
ftst <- fisher.test(t)
pvals <- c(pvals, ftst$p.value)
t <- matrix(c(42,1,55,1), 2, 2, dimnames=list(c("stage I-II", "na"), c("high CD4", "low CD4")))
ftst <- fisher.test(t)
pvals <- c(pvals, ftst$p.value)
p.adjust(pvals, method="BH", n=length(pvals))

count(greaterThan, 'menopause_status')
count(lessThan, 'menopause_status')
t <- matrix(c(13,27,5,16,46,10), 3, 2, dimnames=list(c("pre", "post", "na"), c("high CD4", "low CD4")))
fisher.test(t)

# km DFS
smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$imm
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

summary(smpl.fit, times=1826.25) # 5 year survival
summary(smpl.fit, times=2556.75) # 7 year survival
summary(smpl.fit, times=3652.5, extend=T) # 10 year survival

# cox OS
cox <- coxph(Surv(os_days, os_status)~imm, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

# cox DFS
cox <- coxph(Surv(dfs_days, dfs_status)~imm, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~race, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~pathologic_stage, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~menopause_status, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~imm+pathologic_stage+strata(age_at_initial_pathologic_diagnosis), data=sur)
summary(cox)

# model with all
cox <- coxph(Surv(dfs_days, dfs_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

# remove menopause
cox <- coxph(Surv(dfs_days, dfs_status)~imm+pathologic_stage+race+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

# remove race
cox <- coxph(Surv(dfs_days, dfs_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

# stratified
cox <- coxph(Surv(dfs_days, dfs_status)~imm+strata(pathologic_stage)+race+age_at_initial_pathologic_diagnosis+menopause_status, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~imm+strata(pathologic_stage)+race+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~imm+strata(pathologic_stage), data=sur)
summary(cox)

prop <- cox.zph(cox, transform="km", global=T)

print(prop)
plot(prop)

# km OS
#smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$imm
#sigTest <- survdiff(smpl.surv)
#pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
#smpl.fit <- survfit(smpl.surv)
#pvalue

# cox OS
#cox <- coxph(Surv(os_days, os_status)~imm, data=sur)
#summary(cox)

#cox <- coxph(Surv(os_days, os_status)~age_at_initial_pathologic_diagnosis, data=sur)
#summary(cox)

#----------------------------------
# High CD8/CD4 and Low CD8/CD4
#----------------------------------
clusterGrpFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\clusterGrp.txt"
clusterGrp <- read.delim(clusterGrpFile, header=T, stringsAsFactors=F)
clusterGrp$barcode <- gsub("\\.", "-", clusterGrp$barcode)

sur <- merge(sur, clusterGrp, by.x="bcr_patient_barcode", by.y="barcode", all.x=T)
sur <- sur[!(is.na(sur$clusterGrp)),]

sur <- sur[sur$clusterGrp %in% c("bothNeg", "bothPos"),]
sur$clusterGrp <- factor(sur$clusterGrp, levels=c("bothNeg", "bothPos"), ordered=T)
lowNum <- dim(sur[sur$clusterGrp=="bothNeg",])[1]
highNum <- dim(sur[sur$clusterGrp=="bothPos",])[1]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
#sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

#sur[sur$age_at_initial_pathologic_diagnosis > 30 & sur$age_at_initial_pathologic_diagnosis <= 40,]$age_at_initial_pathologic_diagnosis <- "<=40"
#sur[sur$age_at_initial_pathologic_diagnosis > 40,]$age_at_initial_pathologic_diagnosis <- ">40"
#sur[sur$age_at_initial_pathologic_diagnosis >= 70,]$age_at_initial_pathologic_diagnosis <- ">=70"
#sur$age_at_initial_pathologic_diagnosis <- factor(sur$age_at_initial_pathologic_diagnosis, levels=c("<=40", ">40"), ordered=T)
#count(sur, 'age_at_initial_pathologic_diagnosis')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

greaterThan <- sur[sur$clusterGrp=="bothPos",]
mean(greaterThan$age_at_initial_pathologic_diagnosis)
lessThan <- sur[sur$clusterGrp=="bothNeg",]
mean(lessThan$age_at_initial_pathologic_diagnosis)
wilcox.test(greaterThan$age_at_initial_pathologic_diagnosis, lessThan$age_at_initial_pathologic_diagnosis)

count(greaterThan, 'race')
count(lessThan, 'race')
t <- matrix(c(19,11,3,30,24,4), 3, 2, dimnames=list(c("white", "black", "na"), c("high CD8/CD4", "low CD8/CD4")))
fisher.test(t)

count(greaterThan, 'pathologic_stage')
count(lessThan, 'pathologic_stage')
t <- matrix(c(30,3,0,45,11,2), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("high CD8/CD4", "low CD8/CD4")))
fisher.test(t)

count(greaterThan, 'menopause_status')
count(lessThan, 'menopause_status')
t <- matrix(c(11,19,3,15,34,9), 3, 2, dimnames=list(c("pre", "post", "na"), c("high CD8/CD4", "low CD8/CD4")))
fisher.test(t)

# km OS
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$clusterGrp
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

summary(smpl.fit, times=1826.25) # 5 year survival
summary(smpl.fit, times=2556.75) # 7 year survival
summary(smpl.fit, times=3652.5, extend=T) # 10 year survival

# cox OS
cox <- coxph(Surv(os_days, os_status)~clusterGrp, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~race, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~pathologic_stage, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~menopause_status, data=sur)
summary(cox)

# all
cox <- coxph(Surv(os_days, os_status)~clusterGrp+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~clusterGrp+pathologic_stage+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~clusterGrp+pathologic_stage+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~clusterGrp+pathologic_stage, data=sur)
summary(cox)


# km DFS
#smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$clusterGrp
#sigTest <- survdiff(smpl.surv)
#pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
#smpl.fit <- survfit(smpl.surv)
#pvalue

# cox DFS
cox <- coxph(Surv(dfs_days, dfs_status)~clusterGrp, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~clusterGrp+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~clusterGrp+pathologic_stage+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~clusterGrp+pathologic_stage, data=sur)
summary(cox)

#----------------------------------
# M1 macrophage continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 16
cibersort <- rel[,c(1,cell_idx)]
cellType <- colnames(cibersort)[2]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

# cox os
cox <- coxph(Surv(os_days, os_status)~sum, data=sur)
summary(cox)

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~sum, data=sur)
summary(cox)

#----------------------------------
# M1 macrophages quartiles
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 16
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

quantile(sur$sum)

sur$quartile <- cut(sur$sum, breaks=quantile(sur$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)

q1Num <- dim(sur[sur$quartile=="Q1",])[1]
q2Num <- dim(sur[sur$quartile=="Q2",])[1]
q3Num <- dim(sur[sur$quartile=="Q3",])[1]
q4Num <- dim(sur[sur$quartile=="Q4",])[1]

# km OS
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$quartile
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

#n <- 7
#hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)

# cox os
cox <- coxph(Surv(os_days, os_status)~quartile, data=sur)
summary(cox)

# km DFS
smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$quartile
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
survfit(smpl.surv)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~quartile, data=sur)
summary(cox)

#----------------------------------
# CD8 and CD4 nested coxph categories - not possible to do due to differing numbers between groups
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "cd8"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample")

m <- mean(sur$cd8, na.rm=T)
sd <- sd(sur$cd8, na.rm=T)
times <- 0.25
lowEnd <- round(m - (times * sd), digits=4)
highEnd <- round(m + (times * sd), digits=4)

# lowEnd, highEnd
sur$cd8imm <- NA
sur[sur$cd8 > highEnd & !is.na(sur$cd8),]$cd8imm <- "highcd8"
sur[sur$cd8 < lowEnd & !is.na(sur$cd8),]$cd8imm <- "lowcd8"

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "cd4"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample")

m <- mean(sur$cd4, na.rm=T)
sd <- sd(sur$cd4, na.rm=T)
times <- 0.25
lowEnd <- round(m - (times * sd), digits=4)
highEnd <- round(m + (times * sd), digits=4)

# lowEnd, highEnd
sur$cd4imm <- NA
sur[sur$cd4 > highEnd & !is.na(sur$cd4),]$cd4imm <- "highcd4"
sur[sur$cd4 < lowEnd & !is.na(sur$cd4),]$cd4imm <- "lowcd4"

sur <- sur[!(is.na(sur$cd8imm)),]
sur <- sur[!(is.na(sur$cd4imm)),]

sur$cd8imm <- factor(sur$cd8imm, levels=c("lowcd8", "highcd8"), ordered=T)
sur$cd4imm <- factor(sur$cd4imm, levels=c("lowcd4", "highcd4"), ordered=T)

# covariates
sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)

# km OS
smpl.surv <- Surv(sur$os_days, sur$os_status)~sur$cd8imm
sigTest <- survdiff(smpl.surv)
pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
survfit(smpl.surv)
smpl.fit <- survfit(smpl.surv)
summary(smpl.fit)
pvalue

summary(smpl.fit, times=1826.25) # 5 year survival
summary(smpl.fit, times=2556.75) # 7 year survival
summary(smpl.fit, times=3652.5, extend=T) # 10 year survival

# cox OS
cox <- coxph(Surv(os_days, os_status)~cd8imm, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~race, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~pathologic_stage, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~menopause_status, data=sur)
summary(cox)

# all
coxcd8 <- coxph(Surv(os_days, os_status)~cd8imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(coxcd8)

coxcd8cd4 <- coxph(Surv(os_days, os_status)~cd8imm+cd4imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(coxcd8cd4)

anova(coxcd8,coxcd8cd4)

cox <- coxph(Surv(os_days, os_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status, data=sur)
summary(cox)

cox <- coxph(Surv(os_days, os_status)~imm+pathologic_stage+race+age_at_initial_pathologic_diagnosis, data=sur)
summary(cox)

prop <- cox.zph(cox, transform="km", global=T)

print(prop)
plot(prop)

ggforest(cox, data=sur)

# km DFS
#smpl.surv <- Surv(sur$dfs_days, sur$dfs_status)~sur$imm
#sigTest <- survdiff(smpl.surv)
#pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
#smpl.fit <- survfit(smpl.surv)
#pvalue

# cox DFS
cox <- coxph(Surv(dfs_days, dfs_status)~imm, data=sur)
summary(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~imm+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)


#----------------------------------
# CD8 and cD4 nested coxph continuous
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "cd8"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample")

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "cd4"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample")

#covariates
sur[sur$race=="ASIAN" & !is.na(sur$race),]$race <- NA
sur$race <- factor(sur$race, levels=c(unique(sur$race)))
count(sur, 'race')

sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage I-II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III-IV"
sur[sur$pathologic_stage %in% c("Stage IV"),]$pathologic_stage <- "Stage III-IV"
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I-II", "Stage III-IV"), ordered=T)
count(sur, 'pathologic_stage')

sur[sur$menopause_status %in% c("Indeterminate (neither Pre or Postmenopausal)"),]$menopause_status <- NA
sur[sur$menopause_status %in% c("Peri (6-12 months since last menstrual period)"),]$menopause_status <- NA
sur$menopause_status <- factor(sur$menopause_status, levels=c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                                                              "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"), ordered=T)
count(sur, 'menopause_status')

# cox os
coxcd8 <- coxph(Surv(os_days, os_status)~cd8, data=sur)
summary(coxcd8)
anova(coxcd8)

coxcd8cd4 <- coxph(Surv(os_days, os_status)~cd8+cd4, data=sur)
summary(coxcd8cd4)
anova(coxcd8, coxcd8cd4)

coxcd8 <- coxph(Surv(os_days, os_status)~cd8+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(coxcd8)
anova(coxcd8)

ggforest(coxcd8, data=sur, refLabel="ref")

coxcd8cd4 <- coxph(Surv(os_days, os_status)~cd8+cd4+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(coxcd8cd4)
anova(coxcd8cd4)

anova(coxcd8,coxcd8cd4)

# cox dfs
cox <- coxph(Surv(dfs_days, dfs_status)~cd8, data=sur)
summary(cox)
anova(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~cd4, data=sur)
summary(cox)
anova(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~cd4+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)
anova(cox)

cox <- coxph(Surv(dfs_days, dfs_status)~cd4+cd8+pathologic_stage+age_at_initial_pathologic_diagnosis+menopause_status+race, data=sur)
summary(cox)
anova(cox)

#----------------------------------
# all immune types continuous for purposes of calculating fdr
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

pvalTbl <- data.frame(cellType=character(), survival=character(), splitType=character(), pval=numeric())

for (cell_idx in 2:(dim(rel)[2]-3)){
  
  cibersort <- rel[,c(1,cell_idx)]
  cellType <- colnames(cibersort)[2]
  colnames(cibersort)[2] <- "sum"
  
  surv <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
  
  # analysis of cibersort
  surv <- surv[!(is.na(surv$sum)),]
 
  # cox os
  cox <- coxph(Surv(os_days, os_status)~sum, data=surv)
  summary(cox)
  pvalue <- summary(cox)$coefficients[,5]
  
  comp1 <- c(cellType, "os", "continuous", pvalue)
  pvalTbl[nrow(pvalTbl)+1,] <- comp1
  
  # cox dfs
  cox <- coxph(Surv(dfs_days, dfs_status)~sum, data=surv)
  summary(cox)
  pvalue <- summary(cox)$coefficients[,5]
  
  comp2 <- c(cellType, "dfs", "continuous", pvalue)
  pvalTbl[nrow(pvalTbl)+1,] <- comp2
   
}

cellTblOS <- pvalTbl[pvalTbl$survival=="os",]
cellTblOS$fdr <- p.adjust(cellTblOS$pval, method="fdr")

cellTblDFS <- pvalTbl[pvalTbl$survival=="dfs",]
cellTblDFS$fdr <- p.adjust(cellTblDFS$pval, method="fdr")

#----------------------------------
# all immune types quartiles for purposes of calculating fdr
#----------------------------------
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

pvalTbl <- data.frame(cellType=character(), survival=character(), splitType=character(), pval=numeric())

#cell_idx=2
for (cell_idx in 2:(dim(rel)[2]-3)){
  
  cibersort <- rel[,c(1,cell_idx)]
  cellType <- colnames(cibersort)[2]
  colnames(cibersort)[2] <- "sum"
  
  surv <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
  
  # analysis of cibersort
  surv <- surv[!(is.na(surv$sum)),]
  
  #surv$rank <- rank(surv$sum, ties.method="first")
  #surv$quartile <- cut(surv$rank, breaks=quantile(surv$rank), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)
  #tmp <- surv[order(surv$sum),c(1,75,76,77)]
    
  p <- tryCatch({ # some immune cells have too many 0's such that they can't be split into quartiles
    
    surv$quartile <- cut(surv$sum, breaks=quantile(surv$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)
    
    # km OS
    smpl.surv <- Surv(surv$os_days, surv$os_status)~surv$quartile
    sigTest <- survdiff(smpl.surv)
    pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
    
    }, error = function(e) {
  
    return(NA)
      
    })
  
  comp1 <- c(cellType, "os", "quartile", p)
  pvalTbl[nrow(pvalTbl)+1,] <- comp1

  p <- tryCatch({
    
    surv$quartile <- cut(surv$sum, breaks=quantile(surv$sum), labels=c("Q1","Q2","Q3","Q4"), include.lowest=T)

    # km DFS
    smpl.surv <- Surv(surv$dfs_days, surv$dfs_status)~surv$quartile
    sigTest <- survdiff(smpl.surv)
    pvalue <- 1 - pchisq(sigTest$chisq, length(sigTest$n) - 1)
        
  }, error = function(e) {
    
    return(NA)
    
  })
    
  comp2 <- c(cellType, "dfs", "quartile", p)
  pvalTbl[nrow(pvalTbl)+1,] <- comp2
  
}

cellTblOS <- pvalTbl[pvalTbl$survival=="os",]
cellTblOS$fdr <- p.adjust(cellTblOS$pval, method="fdr")

cellTblDFS <- pvalTbl[pvalTbl$survival=="dfs",]
cellTblDFS$fdr <- p.adjust(cellTblDFS$pval, method="fdr")