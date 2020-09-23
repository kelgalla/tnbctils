#---------------------------------------------------------------------
# FILE     : tnbcType.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-21
# COMMENTS : Take the expression file and place it in a format for 
#            TNBCtype
#            After type is set, see if high vs. low groups are 
#            enriched for a certain subtype
#---------------------------------------------------------------------

library(data.table)
library(plyr)
library(RVAideMemoire)

expFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbcGeneExpSymbol.fpkm.txt"
exp <- read.delim(expFile, header=T, stringsAsFactors=F, na.strings=(""))
exp$mean <- rowMeans(exp[,c(2:133)])

# remove duplicates keeping those with the highest expression mean
exp <- exp[order(exp$hgnc_symbol, -exp$mean),]

dupSym <- unique(exp[duplicated(exp$hgnc_symbol),]$hgnc_symbol)

tmp <- exp[exp$hgnc_symbol %in% dupSym,c(1,134)]

# remove duplicates keeping those with the highest expression mean
exp <- exp[!duplicated(exp$hgnc_symbol),] 

exp <- exp[,c(1:133)]

# tool complains about these being ER+ even though they are not by IHC
exp <- exp[, !(colnames(exp) %in% c("TCGA.BH.A1EW", "TCGA.BH.A0RX", "TCGA.A2.A0ST", "TCGA.A2.A04P", "TCGA.GM.A2DI", 
                               "TCGA.A2.A1G6", "TCGA.E2.A1L7", "TCGA.A2.A3XT"))] 

cm <- as.matrix(exp[,c(2:125)])
rownames(cm) <- exp[,1]

exprs <- cm
sexprs <- scale(exprs)

exp2 <- as.data.frame(sexprs)
setDT(exp2, keep.rownames = "hgnc_symbol")[]

write.csv(exp, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbctype\\tnbcGeneExpSymbol.fpkm.tnbctype.txt",row.names=F)
write.csv(exp2, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbctype\\sexprs\\tnbcGeneExpSymbol.fpkm.tnbctype.sexprs.txt",row.names=F)

#----------------
# after runnning tool, explore subtypes CD8 and CD4

#ttFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbctype\\fd5729e2-5fb9-4a66-9500-1e0b7617e62a_result.csv"
ttFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbctype\\sexprs\\bac566ef-2fce-40d5-a6c0-d8678c8a0cf3_result.csv"

tt <- read.csv(ttFile, header=T, stringsAsFactors=F, na.strings=(""))
colnames(tt)[1] <- "bcr_patient_barcode"
tt$bcr_patient_barcode <- gsub("\\.", "-", tt$bcr_patient_barcode)
tt <- tt[,c(1,2)]
count(tt, 'subtype')

surFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))

rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4
cell_idx <- 8
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
sur <- sur[!(is.na(sur$sum)),]

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

sur <- merge(sur, tt, by.x="bcr_patient_barcode", all.x=T)
sur <- sur[order(sur$imm, sur$subtype),]
sur$subtype <- factor(sur$subtype)

sur[sur$imm==paste0(">", toString(highEnd)),]$subtype
sur[sur$imm==paste0("<", toString(lowEnd)),]$subtype

greaterThan <- sur[sur$imm == CDgreaterThan,]
lessThan <- sur[sur$imm == CDlessThan,]

count(greaterThan, 'subtype')
count(lessThan, 'subtype')

#----------------------------
# CD8 IM and M
t <- matrix(c(3,2,17,3,0,4,6,5, 12,7,1,13,21,4,11,2), 2, 8, byrow=T, 
            dimnames=list(c("high.CD8", "low.CD8"), c("BL1", "BL2", "IM", "LAR", "M", "MSL", "UNS", "na")))
fisher.test(t, workspace = 2e8)
fisher.multcomp(t)

pvals <- c()

t <- matrix(c(3,highNum-3,12,lowNum-12), 2, 2, 
            dimnames=list(c("BL1", "not.BL1"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(2,highNum-2,7,lowNum-7), 2, 2, 
            dimnames=list(c("BL2", "not.BL2"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

# significant
t <- matrix(c(17,highNum-17,1,lowNum-1), 2, 2, 
            dimnames=list(c("IM", "not.IM"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(3,highNum-3,13,lowNum-13), 2, 2, 
            dimnames=list(c("LAR", "not.LAR"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

# significant
t <- matrix(c(0,highNum-0,21,lowNum-21), 2, 2, 
            dimnames=list(c("M", "not.M"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(4,highNum-4,4,lowNum-4), 2, 2, 
            dimnames=list(c("MSL", "not.MSL"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(6,highNum-6,11,lowNum-11), 2, 2, 
            dimnames=list(c("UNS", "not.UNS"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)
fdr <- p.adjust (pvals)

#----------------------------
# CD4 IM and M
t <- matrix(c(7,4,14,4,2,2,9,3, 8,6,4,11,22,7,11,3), 2, 8, byrow=T, 
            dimnames=list(c("high.CD4", "low.CD4"), c("BL1", "BL2", "IM", "LAR", "M", "MSL", "UNS", "na")))
fisher.test(t, workspace = 2e8)
fisher.multcomp(t)

pvals <- c()

t <- matrix(c(7,highNum-7,8,lowNum-8), 2, 2, 
            dimnames=list(c("BL1", "not.BL1"), c("high.CD4", "low.CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(4,highNum-4,6,lowNum-6), 2, 2, 
            dimnames=list(c("BL2", "not.BL2"), c("high.CD4", "low.CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

# significant
t <- matrix(c(14,highNum-14,4,lowNum-4), 2, 2, 
            dimnames=list(c("IM", "not.IM"), c("high.CD4", "low.CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(4,highNum-4,11,lowNum-11), 2, 2, 
            dimnames=list(c("LAR", "not.LAR"), c("high.CD4", "low.CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

# significant
t <- matrix(c(2,highNum-2,22,lowNum-22), 2, 2, 
            dimnames=list(c("M", "not.M"), c("high.CD4", "low.CD4")))
fisher.test(t)

t <- matrix(c(2,highNum-2,7,lowNum-7), 2, 2, 
            dimnames=list(c("MSL", "not.MSL"), c("high.CD4", "low.CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(9,highNum-9,11,lowNum-11), 2, 2, 
            dimnames=list(c("UNS", "not.UNS"), c("high.CD4", "low.CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)
fdr <- p.adjust (pvals)

#-----------
# CD8/CD4 IM and M

#ttFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbctype\\fd5729e2-5fb9-4a66-9500-1e0b7617e62a_result.csv"
ttFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbctype\\sexprs\\bac566ef-2fce-40d5-a6c0-d8678c8a0cf3_result.csv"

tt <- read.csv(ttFile, header=T, stringsAsFactors=F, na.strings=(""))
colnames(tt)[1] <- "bcr_patient_barcode"
tt$bcr_patient_barcode <- gsub("\\.", "-", tt$bcr_patient_barcode)
tt <- tt[,c(1,2)]

surFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))

clusterGrpFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\clusterGrp.txt"
clusterGrp <- read.delim(clusterGrpFile, header=T, stringsAsFactors=F)
clusterGrp$barcode <- gsub("\\.", "-", clusterGrp$barcode)

sur <- merge(sur, clusterGrp, by.x="bcr_patient_barcode", by.y="barcode", all.x=T)
sur <- sur[!(is.na(sur$clusterGrp)),]

sur <- sur[sur$clusterGrp %in% c("bothNeg", "bothPos"),]
sur$clusterGrp <- factor(sur$clusterGrp, levels=c("bothNeg", "bothPos"), ordered=T)
lowNum <- dim(sur[sur$clusterGrp=="bothNeg",])[1]
highNum <- dim(sur[sur$clusterGrp=="bothPos",])[1]

sur <- merge(sur, tt, by.x="bcr_patient_barcode", all.x=T)
sur <- sur[order(sur$clusterGrp, sur$subtype),]
sur$subtype <- factor(sur$subtype)

greaterThan <- sur[sur$clusterGrp=="bothPos",]
lessThan <- sur[sur$clusterGrp=="bothNeg",]

count(greaterThan, 'subtype')
count(lessThan, 'subtype')

t <- matrix(c(3,2,14,3,0,3,4,4, 7,6,1,11,21,3,7,2), 2, 8, byrow=T, 
            dimnames=list(c("high.CD8/CD4", "low.CD8/CD4"), c("BL1", "BL2", "IM", "LAR", "M", "MSL", "UNS", "na")))
fisher.test(t, workspace = 2e8)
fisher.multcomp(t)

pvals <- c()

t <- matrix(c(3,highNum-3,7,lowNum-7), 2, 2, 
            dimnames=list(c("BL1", "not.BL1"), c("high.CD8/CD4", "low.CD8/CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(2,highNum-2,6,lowNum-6), 2, 2, 
            dimnames=list(c("BL2", "not.BL2"), c("high.CD8/CD4", "low.CD8/CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

# significant
t <- matrix(c(14,highNum-14,1,lowNum-1), 2, 2, 
            dimnames=list(c("IM", "not.IM"), c("high.CD8/CD4", "low.CD8/CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(3,highNum-3,11,lowNum-11), 2, 2, 
            dimnames=list(c("LAR", "not.LAR"), c("high.CD8/CD4", "low.CD8/CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

# significant
t <- matrix(c(0,highNum-0,21,lowNum-21), 2, 2, 
            dimnames=list(c("M", "not.M"), c("high.CD8/CD4", "low.CD8/CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(3,highNum-3,3,lowNum-3), 2, 2, 
            dimnames=list(c("MSL", "not.MSL"), c("high.CD8/CD4", "low.CD8/CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(4,highNum-4,7,lowNum-7), 2, 2, 
            dimnames=list(c("UNS", "not.UNS"), c("high.CD8/CD4", "low.CD8/CD4")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)
fdr <- p.adjust (pvals)
