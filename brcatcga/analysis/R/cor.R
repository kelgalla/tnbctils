#---------------------------------------------------------------------
# FILE     : cor.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-17
# COMMENTS : correlation and graphs between tils and cibersort
#---------------------------------------------------------------------

library(ggplot2)
library(ggpubr)

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\"

surFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))

rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

# just lymphoid
#rel$cibertils <- rel$B.cells.naive + rel$B.cells.memory + rel$Plasma.cells + rel$T.cells.CD8 + rel$T.cells.CD4.naive + 
#  rel$T.cells.CD4.memory.resting + rel$T.cells.CD4.memory.activated + rel$T.cells.follicular.helper + 
#  rel$T.cells.regulatory..Tregs. + rel$T.cells.gamma.delta + rel$NK.cells.resting + rel$NK.cells.activated

# all mononuclear cells except PMNs per International Immuno-Oncology Biomarkers Working Group
# left out PMNs and eosinophils which were very low anyway
rel$cibertils <- rel$B.cells.naive + rel$B.cells.memory + rel$Plasma.cells + rel$T.cells.CD8 + rel$T.cells.CD4.naive + 
  rel$T.cells.CD4.memory.resting + rel$T.cells.CD4.memory.activated + rel$T.cells.follicular.helper +
  rel$T.cells.regulatory..Tregs. + rel$T.cells.gamma.delta + rel$NK.cells.resting + rel$NK.cells.activated + rel$Monocytes +
  rel$Macrophages.M0 + rel$Macrophages.M1 + rel$Macrophages.M2 + rel$Dendritic.cells.resting + rel$Dendritic.cells.activated +
  rel$Mast.cells.resting + rel$Mast.cells.activated

# 5 = CD8, 8 = CD4, all = 26, cibertils = 27
cell_idx <- 27
cibersort <- rel[,c(1,cell_idx)]
cellType <- colnames(cibersort)[2]
colnames(cibersort)[2] <- "sum"

surv <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
toCor <- surv[,c("bcr_patient_barcode", "tils", "sum")]
toCor <- toCor[!is.na(toCor$tils) & !is.na(toCor$sum),]

toCor$tilscor <- toCor$tils
toCor[toCor$tilscor %in% c("<1%"),]$tilscor <- 1
toCor[toCor$tilscor %in% c("1-10%"),]$tilscor <- 2
toCor[toCor$tilscor %in% c("10-20%"),]$tilscor <- 3
toCor[toCor$tilscor %in% c("20-30%"),]$tilscor <- 4
toCor[toCor$tilscor %in% c("30-40%"),]$tilscor <- 5
toCor[toCor$tilscor %in% c("40-50%"),]$tilscor <- 6
toCor[toCor$tilscor %in% c("50-60%"),]$tilscor <- 7
#toCor[toCor$tilscor %in% c("60-70%"),]$tilscor <- 8 # no datapoints at this %
toCor[toCor$tilscor %in% c(">70%"),]$tilscor <- 9
toCor$tilscor <- as.numeric(toCor$tilscor)

toCor$tils <- factor(toCor$tils, levels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"), ordered=T)

cor(toCor$tilscor, toCor$sum, method="spearman")
cor.test(toCor$tilscor, toCor$sum, method="spearman")

# cibertils
width <- 2500
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "cor.cibertils.png"), 
    width=width, height=height, res=300)

ggplot(toCor, aes(x=sum, y=tilscor)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(limits=c(0,3), labels=seq(0,3, by=0.5), breaks = seq(0, 3, by = 0.5)) +
  #geom_vline(xintercept=seq(0,3,0.5), colour="white") +
  scale_y_continuous(limits=c(0,10), breaks=c(seq(1, 9, 1)), 
                     labels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"),
                     expand = c(0,0)) +
  xlab("CIBERSORT TILs") +
  ylab("H&E TILs") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  stat_cor(method="spearman", cor.coef.name="'rho'", label.x=0.25, label.y=8.5)

dev.off()

# sample absolute score
width <- 2500
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "cor.png"), 
    width=width, height=height, res=300)

ggplot(toCor, aes(x=sum, y=tilscor)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(limits=c(0,3), labels=seq(0,3, by=0.5), breaks = seq(0, 3, by = 0.5)) +
  #geom_vline(xintercept=seq(0,3,0.5), colour="white") +
  scale_y_continuous(limits=c(0,10), breaks=c(seq(1, 9, 1)), 
                     labels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"),
                     expand = c(0,0)) +
  xlab("CIBERSORT Sample Absolute Score") +
  ylab("TILs") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  stat_cor(method="spearman", cor.coef.name="'rho'", label.x=0.25, label.y=8.5)

dev.off()


# CD8 T cells
width <- 2500
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "cor.T.cells.CD8.png"), 
    width=width, height=height, res=300)

ggplot(toCor, aes(x=sum, y=tilscor)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(limits=c(0,0.7), labels=seq(0,0.7, by=0.1), breaks = seq(0, 0.7, by = 0.1)) +
  #geom_vline(xintercept=seq(0,3,0.5), colour="white") +
  scale_y_continuous(limits=c(0,10), breaks=c(seq(1, 9, 1)), 
                     labels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"),
                     expand = c(0,0)) +
  xlab("CIBERSORT CD8 T Cells") +
  ylab("TILs") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  stat_cor(method="spearman", cor.coef.name="'rho'", label.x=0.35, label.y=8.5)

dev.off()

# CD4 memory activated T cells
width <- 2500
ratio <- 547/600
height <- trunc(width*ratio)
png(paste0(dirSave, "cor.T.cells.CD4.memory.activated.png"), 
    width=width, height=height, res=300)

ggplot(toCor, aes(x=sum, y=tilscor)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(limits=c(0,0.35), labels=seq(0,0.35, by=0.05), breaks = seq(0, 0.35, by = 0.05)) +
  #geom_vline(xintercept=seq(0,3,0.5), colour="white") +
  scale_y_continuous(limits=c(0,10), breaks=c(seq(1, 9, 1)), 
                     labels=c("<1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", ">70%"),
                     expand = c(0,0)) +
  xlab("CIBERSORT CD4 Memory Activated T Cells") +
  ylab("TILs") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  theme(plot.margin=margin(t=5.5, r=35,b=5.5,l=5.5, unit="pt")) +
  stat_cor(method="spearman", cor.coef.name="'rho'", label.x=0.15, label.y=8.5)

dev.off()
