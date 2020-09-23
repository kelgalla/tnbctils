#---------------------------------------------------------------------
# FILE     : mutOfInterest.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2019-07-26
# COMMENTS : Look into the specific mutations in the genes of interest
#---------------------------------------------------------------------

library(plyr)
library(reshape2)

#maf <- read.delim("~/brca/brcatcga/dataReorg/maf/TCGA.BRCA.muse.d565f30b-d3bb-4739-b26c-be8d43c596d4.DR-6.0.protected.maf", header=T, stringsAsFactors=F, comment.char = "#")
maf <- read.delim("F:\\TNBC TILS\\brcatcga\\dataReorg\\maf\\TCGA.BRCA.muse.d565f30b-d3bb-4739-b26c-be8d43c596d4.DR-6.0.protected.maf", header=T, stringsAsFactors=F, comment.char = "#")
#maf <- read.delim("~/brca/brcatcga/dataReorg/maf/TCGA.BRCA.muse.6fcfff20-4993-4789-b27b-69d165130466.DR-6.0.somatic.maf", header=T, stringsAsFactors=F, comment.char = "#")

#tnbc <- read.delim("~/brca/brcatcga/analysis/R/133tnbcCases.txt", header=F, stringsAsFactors=F, col.names=c("tnbc"))
#tnbc <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\133tnbcCases.txt", header=F, stringsAsFactors=F, col.names=c("tnbc"))

#--------------clusterGrp
tnbc <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\clusterGrp.txt", header=T, stringsAsFactors=F)
tnbc <- tnbc[tnbc$clusterGrp=="bothNeg" | tnbc$clusterGrp=="bothPos",]
tnbc$barcode <- gsub("\\.", "-", tnbc$barcode)
colnames(tnbc)[colnames(tnbc)=="barcode"] <- "tnbc"
#--------------

cases <- tnbc$tnbc

#--------------cd8
tnbc <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.txt", header=T, stringsAsFactors=F)
tnbc <- tnbc[!is.na(tnbc$imm),]
colnames(tnbc)[colnames(tnbc)=="barcode"] <- "tnbc"
#--------------

length(tnbc[!(tnbc$tnbc %in% cases),]$tnbc)

cases <- unique(c(cases, tnbc$tnbc))

#--------------cd4
tnbc <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.txt", header=T, stringsAsFactors=F)
tnbc <- tnbc[!is.na(tnbc$imm),]
colnames(tnbc)[colnames(tnbc)=="barcode"] <- "tnbc"
#--------------

length(tnbc[!(tnbc$tnbc %in% cases),]$tnbc)

cases <- unique(c(cases, tnbc$tnbc))

#--------------all cases
tnbc <- data.frame(tnbc=cases)
#---------------

#tnbcMut <- read.delim("~/brca/brcatcga/analysis/bin/tnbcCasesWithMutations.txt", header=F, stringsAsFactors=F, col.names=c("tnbc", "seq"))
tnbcMut <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\bin\\tnbcCasesWithMutations.txt", header=F, stringsAsFactors=F, col.names=c("tnbc", "seq"))

tnbcMut[tnbcMut$seq=="mutation found",]$seq <- 1
tnbcMut[tnbcMut$seq=="no mutations",]$seq <- 0

tnbc <- merge(tnbc,tnbcMut, by="tnbc")

dat <- data.frame(Hugo_Symbol=maf$Hugo_Symbol, Entrez_Gene_Id=maf$Entrez_Gene_Id, Tumor_Sample_Barcode=maf$Tumor_Sample_Barcode, 
                  Variant_Classification=maf$Variant_Classification, IMPACT=maf$IMPACT, Mutation_Status=maf$Mutation_Status,
                  HGVSc=maf$HGVSc, HGVSp=maf$HGVSp, HGVSp_Short=maf$HGVSp_Short, stringsAsFactors=F)

#temp <- dat[dat$Variant_Classification=="Missense_Mutation",] # keep variants that have disruptive impact on the protein
dat <- dat[dat$IMPACT %in% c("HIGH","MODERATE"),] # keep variants that have disruptive impact on the protein

dat$barcode <- gsub("(TCGA-[a-zA-Z0-9]{2}-[a-zA-Z0-9]{4}).*", "\\1", dat$Tumor_Sample_Barcode, perl=T)

dat$mutation <- 1

dat <- dat[dat$barcode %in% tnbc$tnbc,] # be sure to only keep our subset of interest

datUniq <- ddply(dat, .(Hugo_Symbol, barcode), summarize, Variant_Classification=paste(Variant_Classification, collapse=" | "),
                 HGVSc=paste(HGVSc, collapse=" | "), HGVSp=paste(HGVSp, collapse=" | "), HGVSp_Short=paste(HGVSp_Short, collapse=" | "),
                 mutation = 1) 

#---------clusterGrp
goi <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.clusterGrp.mutSigDiff.txt", header=T, stringsAsFactors=F)
#---------

genes <- goi$Hugo_Symbol

#---------cd8
goi <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
#---------

length(goi[!(goi$Hugo_Symbol %in% genes),]$Hugo_Symbol)

genes <- unique(c(genes, goi$Hugo_Symbol))

#---------cd4
goi <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
#---------

length(goi[!(goi$Hugo_Symbol %in% genes),]$Hugo_Symbol)

genes <- unique(c(genes, goi$Hugo_Symbol))

#---------all
goi <- data.frame(Hugo_Symbol=genes)
#---------

datInterest <- datUniq[datUniq$Hugo_Symbol %in% goi$Hugo_Symbol,]

colnames(datInterest)[colnames(datInterest)=="barcode"] <- "Sample_ID"
colnames(datInterest)[colnames(datInterest)=="HGVSp_Short"] <- "Protein_Change"

#------------clusterGrp
for (i in unique(datInterest$Hugo_Symbol)){
  
write.table(datInterest[datInterest$Hugo_Symbol==i, c("Hugo_Symbol", "Sample_ID", "Protein_Change")], 
            paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mutations\\clusterGrp.", i, ".txt"), sep="\t", row.names=F, quote=F)
}
#-------------

#------------cd8
for (i in unique(datInterest$Hugo_Symbol)){
  
  write.table(datInterest[datInterest$Hugo_Symbol==i, c("Hugo_Symbol", "Sample_ID", "Protein_Change")], 
              paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mutations\\cd8.", i, ".txt"), sep="\t", row.names=F, quote=F)
} 
#-------------

#-------------cd4
for (i in unique(datInterest$Hugo_Symbol)){
  
  write.table(datInterest[datInterest$Hugo_Symbol==i, c("Hugo_Symbol", "Sample_ID", "Protein_Change")], 
              paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mutations\\cd4.", i, ".txt"), sep="\t", row.names=F, quote=F)
} 
#-------------

#-------------cd4
for (i in unique(datInterest$Hugo_Symbol)){
  
  write.table(datInterest[datInterest$Hugo_Symbol==i, c("Hugo_Symbol", "Sample_ID", "Protein_Change")], 
              paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mutations\\all.", i, ".txt"), sep="\t", row.names=F, quote=F)
} 
#-------------