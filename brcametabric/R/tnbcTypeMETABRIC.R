#---------------------------------------------------------------------
# FILE     : tnbcTypeMETABRIC.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-21
# COMMENTS : Take the expression file and place it in a format for 
#            TNBCtype
#            After type is set, see if high vs. low groups are 
#            enriched for a certain subtype
#---------------------------------------------------------------------

library(dplyr)
library(RVAideMemoire)

# going to use the intensity values instead of log2 intensity as that is what was used for CIBERSORT
expFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\expForCIBERSORT.txt"
exp <- read.delim(expFile, header=T, stringsAsFactors=F, na.strings=(""))

# tool complains about these being ER+ even though they are not
exp <- exp[, !(colnames(exp) %in% c("MB.3297", "MB.7057", "MB.2993", "MB.7087", "MB.5008", "MB.6052", "MB.4303", "MB.7031", "MB.5155"))]

cm <- as.matrix(exp[,c(2:188)])
rownames(cm) <- exp[,1]

exprs <- cm
sexprs <- scale(exprs)

exp2 <- as.data.frame(sexprs)
setDT(exp2, keep.rownames = "Hugo_Symbol")[]

write.csv(exp, "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\expForCIBERSORT.tnbctype.txt",row.names=F)
write.csv(exp2, "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\expForCIBERSORT.tnbctype.sexprs.txt",row.names=F)

#-------------
# after runnning tool, explore subtypes

#ttFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\tnbctype\\1b17c667-4503-4db3-a057-ff4dcc183b5b_result.csv"
ttFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\tnbctype\\sexprs\\1728e906-266b-473f-b314-97ea92023e10_result.csv"

tt <- read.csv(ttFile, header=T, stringsAsFactors=F, na.strings=(""))
colnames(tt)[1] <- "PATIENT_ID"
tt$PATIENT_ID <- gsub("\\.", "-", tt$PATIENT_ID)
tt <- tt[,c(1,2)]
count(tt, "subtype")

surFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcSurv.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
sur$OS_STATUS2 <- ifelse(sur$OS_STATUS == "0:LIVING", 0, 1)

rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# 5 = CD8, 8 = CD4, 11 = gd
cell_idx <- 11
cibersort <- rel[,c(1,cell_idx)]
colnames(cibersort)[2] <- "sum"

sur <- merge(sur, cibersort, by.x="PATIENT_ID", by.y="Input.Sample", all.x=T)
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

sur <- merge(sur, tt, by.x="PATIENT_ID", all.x=T)
sur <- sur[order(sur$imm, sur$subtype),]

sur[sur$imm==paste0(">", toString(highEnd)),]$subtype
sur[sur$imm==paste0("<", toString(lowEnd)),]$subtype

greaterThan <- sur[sur$imm == CDgreaterThan,]
lessThan <- sur[sur$imm == CDlessThan,]

count(greaterThan, 'subtype')
count(lessThan, 'subtype')

#----------------------------
# CD8 T cells
t <- matrix(c(7,4,26,5,6,4,10,3, 15,8,9,12,31,3,18,5), 2, 8, byrow=T,
            dimnames=list(c("high.CD8", "low.CD8"), c("BL1", "BL2", "IM", "LAR", "M", "MSL", "UNS", "na")))
fisher.test(t, workspace = 2e8)
fisher.multcomp(t)

pvals <- c()

t <- matrix(c(7,highNum-7,15,lowNum-15), 2, 2, 
            dimnames=list(c("BL1", "not.BL1"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(4,highNum-4,8,lowNum-8), 2, 2, 
            dimnames=list(c("BL2", "not.BL2"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(26,highNum-26,9,lowNum-9), 2, 2, 
            dimnames=list(c("IM", "not.IM"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(5,highNum-5,12,lowNum-12), 2, 2, 
            dimnames=list(c("LAR", "not.LAR"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(6,highNum-6,31,lowNum-31), 2, 2, 
            dimnames=list(c("M", "not.M"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(4,highNum-4,3,lowNum-3), 2, 2, 
            dimnames=list(c("MSL", "not.MSL"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(10,highNum-10,18,lowNum-18), 2, 2, 
            dimnames=list(c("UNS", "not.UNS"), c("high.CD8", "low.CD8")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)
fdr <- p.adjust(pvals)

#----------------------------
# gd T cells
t <- matrix(c(10,4,17,8,6,6,12,3, 15,9,19,9,27,4,15,5), 2, 8, byrow=T,
            dimnames=list(c("high.gdT", "low.gdT"), c("BL1", "BL2", "IM", "LAR", "M", "MSL", "UNS", "na")))
fisher.test(t, workspace = 2e8)
fisher.multcomp(t)

pvals <- c()

t <- matrix(c(10,highNum-10,13,lowNum-13), 2, 2, 
            dimnames=list(c("BL1", "not.BL1"), c("high.gdT", "low.gdT")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(7,highNum-7,6,lowNum-6), 2, 2, 
            dimnames=list(c("BL2", "not.BL2"), c("high.gdT", "low.gdT")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(17,highNum-17,17,lowNum-17), 2, 2, 
            dimnames=list(c("IM", "not.IM"), c("high.gdT", "low.gdT")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(7,highNum-7,8,lowNum-8), 2, 2, 
            dimnames=list(c("LAR", "not.LAR"), c("high.gdT", "low.gdT")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

# significant
t <- matrix(c(6,highNum-6,27,lowNum-27), 2, 2, 
            dimnames=list(c("M", "not.M"), c("high.gdT", "low.gdT")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(6,highNum-6,5,lowNum-5), 2, 2, 
            dimnames=list(c("MSL", "not.MSL"), c("high.gdT", "low.gdT")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)

t <- matrix(c(10,highNum-10,22,lowNum-22), 2, 2, 
            dimnames=list(c("UNS", "not.UNS"), c("high.gdT", "low.gdT")))
fisher.test(t)
p <- fisher.test(t)$p.value

pvals <- append(pvals, p)
fdr <- p.adjust(pvals)

