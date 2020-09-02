#---------------------------------------------------------------------
# FILE     : orderedMutHeatmap.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2018-02-28
# COMMENTS : this script will read in data from a bunch of different 
#            files for the purpose of generating a subtype specific mutation heatmap
#---------------------------------------------------------------------

#dir <- "~/brca/brcatcga/analysis/R"
dir <- "F:\\TNBC TILS\\brcatcga\\analysis\\R"

# import gene list of interest to plot

#------- CD8 T cells
#mutSigDiff <- read.delim("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD8.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
mutSigDiff <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
#mutSigDiff <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.mediansd.mutSigDiff.txt", header=T, stringsAsFactors=F)

#------- CD4 memory activated T cells
#mutSigDiff <- read.delim("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
mutSigDiff <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
#mutSigDiff <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.mediansd.mutSigDiff.txt", header=T, stringsAsFactors=F)

#------- CD8 CD4
#mutSigDiff <- read.delim("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
mutSigDiff <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.clusterGrp.mutSigDiff.txt", header=T, stringsAsFactors=F)

genes <- mutSigDiff$Hugo_Symbol

# import immune status of sample for sorting
#------- CD8 T cells
#pts <- read.delim(paste0(dir, "/os.abs.T.cells.CD8.0.25sd.txt"), header=T, stringsAsFactors=F)
pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD8.0.25sd.txt"), header=T, stringsAsFactors=F)
#pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD8.mediansd.txt"), header=T, stringsAsFactors=F)
#------- CD4 memory activated T cells
#pts <- read.delim(paste0(dir, "/os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), header=T, stringsAsFactors=F)
pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), header=T, stringsAsFactors=F)
#pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD4.memory.activated.mediansd.txt"), header=T, stringsAsFactors=F)

#------- CD8 CD4
#pts <- read.delim(paste0(dir, "/os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), header=T, stringsAsFactors=F)
pts <- read.delim(paste0(dir, "\\clusterGrp.txt"), header=T, stringsAsFactors=F)
pts <- pts[pts$clusterGrp %in% c("bothPos", "bothNeg"),]
pts[pts$clusterGrp=="bothPos",]$clusterGrp <- "high"
pts[pts$clusterGrp=="bothNeg",]$clusterGrp <- "low"
colnames(pts)[colnames(pts)=="clusterGrp"] <- "imm"

pts$barcode <- gsub("-", "\\.", pts$barcode)

# import mutation count of sample for sorting
#mutCnt <- read.delim("~/brca/brcatcga/analysis/R/mutCnt.txt", header=T, stringsAsFactors=F)
mutCnt <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mutCnt.txt", header=T, stringsAsFactors=F)

# combine immune and mutation data
mutSum <- merge(mutCnt, pts, by.x="patient", by.y="barcode")
mutSum <- mutSum[!is.na(mutSum$imm),] # remove imm NAs
mutSum <- mutSum[!is.na(mutSum$mutCnt),] # remove mutCnt NAs
mutSum <- mutSum[order(mutSum$imm, mutSum$mutCnt),]

# get non-summarized mutation data for graphing
#datNew <- read.delim("~/brca/brcatcga/analysis/R/mut.txt", header=T, stringsAsFactors=F)
datNew <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mut.txt", header=T, stringsAsFactors=F)
datNew$barcode <- gsub("-", "\\.", datNew$barcode)
datp <- datNew

datp <- datp[datp$Hugo_Symbol %in% genes,] # get subset of genes to graph
datp$Hugo_Symbol <- factor(datp$Hugo_Symbol, rev(genes))

datp <- datp[datp$barcode %in% mutSum$patient,] # get subset of samples to graph
datp$barcode <- factor(datp$barcode, levels=mutSum$patient) # put in desired order

# dummy column for mutation data so I can color on the subtypes, 2 = high mutant, 3 = low mutant
datp$modMut <- ifelse(datp$mutation==0, 0, ifelse(datp$barcode %in% mutSum[mutSum$imm=="high",]$patient, 2 , 3))
datp$mutation <- factor(datp$mutation)
datp$modMut <- factor(datp$modMut)

width <- 2500
#ratio <- 461/572
ratio <- 547/600
height <- trunc(width*ratio)
#---- CD8 T Cells
#png(paste0("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD8.0.25sd.mutHeatmap2.png"), width=width, height=height)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.mutHeatmap2.fixline.png"), width=width, height=height, res=300)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.mediansd.mutHeatmap2.png"), width=width, height=height)

#---- CD4 memory activated T Cells
#png(paste0("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutHeatmap2.png"), width=width, height=height)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutHeatmap2.fixline.png"), width=width, height=height, res=300)
#png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.mediansd.mutHeatmap2.png"), width=width, height=height)

#---- CD8 CD4
#png(paste0("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutHeatmap2.png"), width=width, height=height)
png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.clusterGrp.mutHeatmap2.fixline.png"), width=width, height=height, res=300)

n1 <- length(unique(datp$Hugo_Symbol))
n2 <- length(unique(datp$barcode))

# sig genes
ggplot(datp) + geom_tile(aes(x=barcode, y=Hugo_Symbol, fill=modMut)) +
  #geom_tile(height=0.97) +
  #scale_fill_manual(values = c("#FFFFFF","#323D8D"), na.value="#BEBEBE") +
  #scale_x_discrete(expand=c(0, 1, 0, 1)) +
  #scale_y_discrete(expand=c(0, 0, 0, 0)) +
  scale_x_discrete(expand = expand_scale(add=0)) +
  scale_y_discrete(expand = expand_scale(add=0)) +
  guides(fill=F) +
  ylab("Genes") +
  xlab("Sample") +
  #scale_fill_manual(values=c("#FFFFFF", "#00BFC4", "#F8766D"),na.value="#BEBEBE") +
  scale_fill_manual(values=c("#FFFFFF", "#EB00FFFF", "#FF9900FF"), na.value="#BEBEBE") + # cluster
  #ggtitle("Top 34 TNBC mutations") +
  #theme(text = element_text(size=14, colour="black")) +
  #theme(axis.text = element_text(size=12, colour="black")) +
  #theme(axis.text.x = element_text(size=8, vjust=0.5, hjust=1, angle=90, colour="black")) +
  theme(axis.text.y = element_text(size=20, colour="black", face="bold")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(panel.grid.minor.x=element_blank()) +
  theme(panel.border=element_rect(colour="black", fill=NA, size=2)) +
  theme(panel.grid.major.x=element_blank()) +
  theme(panel.grid.major.y=element_blank()) +
  theme(panel.grid.minor.y=element_blank()) +
  theme(panel.background=element_rect(fill="#000000")) +
  geom_line(data = data.frame(x = c(0, n2) + 0.5, y = rep(2:n1, each = 2) - 0.5), aes(x = x, y = y, group = y))

dev.off()

#ggsave("os.abs.T.cells.CD8.0.25sd.mutHeatmap2.fixnoline.png", plot=plot, path=paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R"),
#       width=8.333, height=7.597, units="in", dpi=300)
