#---------------------------------------------------------------------
# FILE     : mutFreq.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2018-02-19
# COMMENTS : read the mutTable.txt file and calculate the mutation 
#            frequency for different subsets of tnbc cases and save 
#            it in a specific version of mutTableNew.txt (count in 
#            last row)
#            
#            save the general mutation count for all the tnbc cases in 
#            mutCnt.txt
#---------------------------------------------------------------------

library(ggplot2)

#dir <- "~/brca/brcatcga/analysis/R"
dir <- "F:\\TNBC TILS\\brcatcga\\analysis\\R"

#mutTable <- read.delim(paste0(dir, "/mutTable.txt"), header=T, stringsAsFactors=F)
mutTable <- read.delim(paste0(dir, "\\mutTable.txt"), header=T, stringsAsFactors=F)
gns <- dim(mutTable)[1]

#------------
# CD8 T Cells
#pts <- read.delim(paste0(dir, "/os.abs.T.cells.CD8.0.25sd.txt"), header=T, stringsAsFactors=F)
#pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD8.0.25sd.txt"), header=T, stringsAsFactors=F)
#pts$barcode <- gsub("-", "\\.", pts$barcode)

# median
#pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD8.mediansd.txt"), header=T, stringsAsFactors=F)
#pts$barcode <- gsub("-", "\\.", pts$barcode)
#-----------

#------------
# CD4 Memory Activated T Cells
#pts <- read.delim(paste0(dir, "/os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), header=T, stringsAsFactors=F)
pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), header=T, stringsAsFactors=F)
pts$barcode <- gsub("-", "\\.", pts$barcode)

#pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD4.memory.activated.mediansd.txt"), header=T, stringsAsFactors=F)
#pts$barcode <- gsub("-", "\\.", pts$barcode)
#------------

#------------
# CD8+CD4
#pts <- read.delim(paste0(dir, "/os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), header=T, stringsAsFactors=F)
#pts$barcode <- gsub("-", "\\.", pts$barcode)

pts <- read.delim(paste0(dir, "\\clusterGrp.txt"), header=T, stringsAsFactors=F)
pts$barcode <- gsub("-", "\\.", pts$barcode)
pts <- pts[pts$clusterGrp %in% c("bothPos", "bothNeg"),]
pts[pts$clusterGrp=="bothPos",]$clusterGrp <- "high"
pts[pts$clusterGrp=="bothNeg",]$clusterGrp <- "low"
colnames(pts)[colnames(pts)=="clusterGrp"] <- "imm"
#------------

high <- pts[pts$imm=="high" & !is.na(pts$imm),c("barcode")]
highMut <- mutTable[,colnames(mutTable) %in% high]

low <- pts[pts$imm=="low" & !is.na(pts$imm),c("barcode")]
lowMut <- mutTable[,colnames(mutTable) %in% low]

mutTable$highNum <-  apply(highMut, 1, function (x) {int <- as.integer(x) 
                                                sum(as.integer(int), na.rm=T)})
mutTable$highCases <-  apply(highMut, 1, function (x) {int <- as.integer(x)
                                                  length(int[!is.na(int)])})
mutTable$highFreq <-  apply(highMut, 1, function (x) {int <- as.integer(x)
                                                 (sum(int, na.rm=T)/length(int[!is.na(int)]))*100})

mutTable$lowNum <-  apply(lowMut, 1, function (x) {int <- as.integer(x) 
                                                  sum(as.integer(int), na.rm=T)})
mutTable$lowCases <-  apply(lowMut, 1, function (x) {int <- as.integer(x)
                                                    length(int[!is.na(int)])})
mutTable$lowFreq <-  apply(lowMut, 1, function (x) {int <- as.integer(x)
                                                   (sum(int, na.rm=T)/length(int[!is.na(int)]))*100})

mutSum <- apply(mutTable[,2:134], 2, function (x) sum(as.numeric(x)))
mutCnt <- data.frame(patient=colnames(mutTable[,2:134]), mutCnt=mutSum, stringsAsFactors=F)
median(mutCnt$mutCnt, na.rm=T)
mean(mutCnt$mutCnt, na.rm=T)

#write.table(mutCnt, "~/brca/brcatcga/analysis/R/mutCnt.txt", sep="\t", row.names=F, quote=F)

mutTotal <- c("Total", mutSum, NA, NA, NA, NA, NA, NA, NA, NA, NA)
mutTableNew <- rbind(mutTable, mutTotal)

#------
# CD8 T Cells
#write.table(mutTableNew, "~/brca/brcatcga/analysis/R/os.abs.T.cells.CD8.0.25sd.mutTable.txt", sep="\t", row.names=F, quote=F)
#write.table(mutTableNew, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.mediansd.mutTable.txt", sep="\t", row.names=F, quote=F)
#---------

#---------
# CD4 memory activated T cells
#write.table(mutTableNew, "~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutTable.txt", sep="\t", row.names=F, quote=F)
#write.table(mutTableNew, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.mediansd.mutTable.txt", sep="\t", row.names=F, quote=F)
#---------

#---------
# CD8 CD4
#write.table(mutTableNew, "~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutTable.txt", sep="\t", row.names=F, quote=F)
#write.table(mutTableNew, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.clusterGrp.mutTable.txt", sep="\t", row.names=F, quote=F)
#---------

interest <- mutTableNew[,c(1,135:143)]

ok <- tail(mutTableNew)

# graph the frequency of gene mutations overall for tnbc
mutTableFreq <- mutTable[,c("Hugo_Symbol", "freq")]
mutTableFreq$Hugo_Symbol <- factor(mutTableFreq$Hugo_Symbol, levels=mutTableFreq$Hugo_Symbol)

width <- 2940
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
#png(paste0("~/brca/brcatcga/analysis/R/top34mut.png"), width=width, height=height, res=300)
png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\top34mutSize.png"), width=width, height=height, res=300)

ggplot(data=mutTableFreq[1:34,], aes(x=Hugo_Symbol, y=freq, fill=Hugo_Symbol)) + 
  geom_bar(stat="identity", colour="black") +
  guides(fill=F) +
  xlab("Gene") +
  #ylab("Mutation Frequency") +
  #ggtitle("Frequency of Gene Mutations") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=14, colour="black")) +
  theme(axis.text.x = element_text(size=20, angle=90, vjust=0.25, hjust=1)) +
  theme(axis.text.y = element_text(size=20, angle=90, hjust=0.5)) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  #theme(legend.text = element_text(size=14, colour="black")) +  
  #theme(legend.title = element_text(size=16, colour="black")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  scale_y_continuous(limits=c(0, 80)) +
  geom_text(aes(label=paste0(format(freq,digits=2), "%"), angle=90, vjust=0.4, hjust=-0.1, size=20)) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank())

dev.off()

# look for differences in mutation count between the groups, and graph the mutation count in general for all cases
mutSum <- merge(mutCnt, pts, by.x="patient", by.y="barcode")

mutSumDF <- mutSum[!is.na(mutSum$imm),] # remove imm NAs
mutSumDF <- mutSumDF[!is.na(mutSumDF$mutCnt),] # remove mutCnt NAs

mutSumDF$imm <- factor(mutSumDF$imm, levels=c("high", "low"))
mutSumDF <- mutSumDF[order(mutSumDF$imm, mutSumDF$mutCnt),]
mutSumDF$patient <- factor(mutSumDF$patient, levels=mutSumDF$patient)

median(mutSumDF[mutSumDF$imm=="high",c("mutCnt")])
median(mutSumDF[mutSumDF$imm=="low",c("mutCnt")])
mean(mutSumDF[mutSumDF$imm=="high",c("mutCnt")])
mean(mutSumDF[mutSumDF$imm=="low",c("mutCnt")])

#summary(aov(mutSumDF$mutCnt ~ mutSumDF$imm))
#TukeyHSD(aov(mutSumDF$mutCnt ~ mutSumDF$imm))
t.test(mutSumDF$mutCnt ~ mutSumDF$imm)
wilcox.test(mutSumDF$mutCnt ~ mutSumDF$imm)
#kruskal.test(mutSumDF$mutCnt ~ mutSumDF$imm)
#pairwise.wilcox.test(mutSumDF$mutCnt, mutSumDF$imm, p.adj="bonferroni")

# graph mutation count for each of the samples in the different groups

width <- 2940
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
#png(paste0("~/brca/brcatcga/analysis/R/mutCnt.png"), width=width, height=height)
png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mutCntSize.png"), width=width, height=height, res=300)

# plot without high vs. low
#ggplot(data=mutSumDF, aes(x=patient, y=mutCnt, fill=imm)) +
ggplot(data=mutSum, aes(x=patient, y=mutCnt)) +
#geom_bar(stat="identity", colour="black") +
  geom_bar(stat="identity", fill="#323D8D", colour="black") +
  #geom_bar(stat="identity", colour="black", fill=c("#D39200", "#F8766D", "#00BA38", "#619CFF"))
  #geom_bar(stat="identity", colour="black", fill=c(rep("#D39200", 3), rep("#F8766D", 40), 
  #      rep("#00BA38", 38), rep("#619CFF", 19))) +
  #guides(fill=F) +
  #scale_fill_manual(values=c("#00BFC4", "#F8766D"), 
   #                 name="CD8 T Cells",
  #                  breaks=c("high", "low"),
  #                  labels=c("High", "Low")) +
  xlab("Sample") +
  ylab("Number of mutations per sample") +
  #ggtitle("Number of Mutations") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.text.x = element_text(size=8, angle=-90, vjust=0.25, hjust=0)) +
  theme(axis.text.y = element_text(size=20, colour="black", face="bold")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(legend.text = element_text(size=12, colour="black")) +  
  theme(legend.title = element_text(size=14, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(legend.key = element_rect(color="black", size=1.25)) +
  #geom_text(aes(label=mutCnt, angle=90, vjust=0.4, hjust=-0.1), size=3.25) +
  theme(axis.text.x=element_blank())


  #scale_y_continuous(limits=c(0, 150), breaks=c(0, 25, 50, 75, 100, 125, 150))

dev.off()

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)

#---- CD8 T cells
#cellType <- "CD8 T Cells"
#png(paste0("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD8.0.25sd.mutCnt.png"), width=width, height=height)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.mutCnt.png"), width=width, height=height, res=300)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.mediansd.mutCnt.png"), width=width, height=height, res=300)

#---- CD4 memory activated T cells
cellType <- "CD4 Memory\nActivated T Cells"
#png(paste0("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutCnt.png"), width=width, height=height)
png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutCnt.png"), width=width, height=height, res=300)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.mediansd.mutCnt.png"), width=width, height=height, res=300)

#---- CD8 CD4
#cellType <- "CD8 T Cells +\nCD4 Memory\nActivated T Cells"
#png(paste0("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutCnt.png"), width=width, height=height)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.clusterGrp.mutCnt.png"), width=width, height=height, res=300)

# plot with legend
ggplot(data=mutSumDF, aes(x=patient, y=mutCnt, fill=imm)) +
  #geom_bar(stat="identity", colour="black") +
  geom_bar(stat="identity", colour="black") +
  #geom_bar(stat="identity", colour="black", fill=c("#D39200", "#F8766D", "#00BA38", "#619CFF"))
  #geom_bar(stat="identity", colour="black", fill=c(rep("#D39200", 3), rep("#F8766D", 40), 
  #      rep("#00BA38", 38), rep("#619CFF", 19))) +
  #guides(fill=F) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D"), name=cellType, breaks=c("high", "low"), labels=c("High", "Low")) +
  #scale_fill_manual(values=c("#EB00FFFF", "#FF9900FF"), name=cellType, breaks=c("high", "low"), labels=c("High", "Low")) + # cluster
  xlab("Sample") +
  ylab("Number of mutations per sample") +
  #ggtitle("Number of Mutations") +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=20, colour="black")) +
  theme(axis.text.x = element_text(size=8, angle=-90, vjust=0.25, hjust=0)) +
  theme(axis.text.y = element_text(size=20, colour="black", face="bold")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(legend.text = element_text(size=18, colour="black")) +  
  theme(legend.title = element_text(size=20, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(legend.key = element_rect(color="black", size=1.25)) +
  #geom_text(aes(label=mutCnt, angle=90, vjust=0.4, hjust=-0.1), size=3.25) +
  theme(axis.text.x=element_blank()) +
  theme(legend.position = c(0.7, 0.85))

# theme.size = (14/5) * geom.text.size

dev.off()
