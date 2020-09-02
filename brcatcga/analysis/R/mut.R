#---------------------------------------------------------------------
# FILE     : mut.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2018-02-19
# COMMENTS : read the maf file, and summarize each gene as being 
#            mutated or not for the tnbc cases, save the data in 
#            mut.txt and mutTable.txt, also graph the data for the 
#            top 34 genes in a mutation heatmap
#---------------------------------------------------------------------

library(plyr)
library(reshape2)

#maf <- read.delim("~/brca/brcatcga/dataReorg/maf/TCGA.BRCA.muse.d565f30b-d3bb-4739-b26c-be8d43c596d4.DR-6.0.protected.maf", header=T, stringsAsFactors=F, comment.char = "#")
maf <- read.delim("F:\\TNBC TILS\\brcatcga\\dataReorg\\maf\\TCGA.BRCA.muse.d565f30b-d3bb-4739-b26c-be8d43c596d4.DR-6.0.protected.maf", header=T, stringsAsFactors=F, comment.char = "#")
#maf <- read.delim("~/brca/brcatcga/dataReorg/maf/TCGA.BRCA.muse.6fcfff20-4993-4789-b27b-69d165130466.DR-6.0.somatic.maf", header=T, stringsAsFactors=F, comment.char = "#")

#tnbc <- read.delim("~/brca/brcatcga/analysis/R/133tnbcCases.txt", header=F, stringsAsFactors=F, col.names=c("tnbc"))
tnbc <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\133tnbcCases.txt", header=F, stringsAsFactors=F, col.names=c("tnbc"))
#tnbcMut <- read.delim("~/brca/brcatcga/analysis/bin/tnbcCasesWithMutations.txt", header=F, stringsAsFactors=F, col.names=c("tnbc", "seq"))
tnbcMut <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\bin\\tnbcCasesWithMutations.txt", header=F, stringsAsFactors=F, col.names=c("tnbc", "seq"))
tnbcMut[tnbcMut$seq=="mutation found",]$seq <- 1
tnbcMut[tnbcMut$seq=="no mutations",]$seq <- 0

tnbc <- merge(tnbc,tnbcMut, by="tnbc")

dat <- data.frame(Hugo_Symbol=maf$Hugo_Symbol, Entrez_Gene_Id=maf$Entrez_Gene_Id, Tumor_Sample_Barcode=maf$Tumor_Sample_Barcode, 
                  Variant_Classification=maf$Variant_Classification, IMPACT=maf$IMPACT, Mutation_Status=maf$Mutation_Status, stringsAsFactors=F)

#temp <- dat[dat$Variant_Classification=="Missense_Mutation",] # keep variants that have disruptive impact on the protein
dat <- dat[dat$IMPACT %in% c("HIGH","MODERATE"),] # keep variants that have disruptive impact on the protein

dat$barcode <- gsub("(TCGA-[a-zA-Z0-9]{2}-[a-zA-Z0-9]{4}).*", "\\1", dat$Tumor_Sample_Barcode, perl=T)

dat$mutation <- 1

dat <- dat[dat$barcode %in% tnbc$tnbc,] # be sure to only keep our subset of interest
datUniq <- ddply(dat, .(Hugo_Symbol, barcode), summarize, Variant_Classification=paste(Variant_Classification, collapse=" | "), mutation = 1) 

datSum <- ddply(datUniq, .(Hugo_Symbol), summarize, mutNum = length(mutation)) 
datSum <- datSum[rev(order(datSum$mutNum)),]
#genes <- datSum[1:8453,c("Hugo_Symbol")]
genes <- datSum[1:34,c("Hugo_Symbol")]
genes <- as.character(rev(genes))
#genes <- data.frame(genes=genes, gorder=1:length(genes))

datNew <- datUniq

# takes a really long time

allgenes <- datSum[,c("Hugo_Symbol")]
allgenes <- as.character(rev(allgenes))
for (n in allgenes) {
  #n <- "TP53"
  sub <- datNew[datNew$Hugo_Symbol==n,]
  nomut <- tnbc[!(tnbc$tnbc %in% sub$barcode),]
  nomutseq <- nomut[nomut$seq==1,]
  nomutseqDF <- data.frame(Hugo_Symbol=n, barcode=nomutseq$tnbc, Variant_Classification="WT", mutation=0)
  noseq <- nomut[nomut$seq==0,]
  noseqDF <- data.frame(Hugo_Symbol=n, barcode=noseq$tnbc, Variant_Classification=NA, mutation=NA)
  datNew <- rbind(datNew, nomutseqDF, noseqDF)
  #datNew <- rbind(datNew, nomutseqDF)
}
tmp <- datNew[datNew$Hugo_Symbol=="TP53",]
write.table(datNew, "~/brca/brcatcga/analysis/R/mut.txt", sep="\t", row.names=F, quote=F)

#datNew <- read.delim("~/brca/brcatcga/analysis/R/mut.txt", header=T, stringsAsFactors=F)
datNew <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\mut.txt", header=T, stringsAsFactors=F)
datp <- datNew

datp <- datp[datp$Hugo_Symbol %in% genes,] # get subset to graph
datp$barcode <- factor(datp$barcode, tnbc$tnbc)
datp$Hugo_Symbol <- factor(datp$Hugo_Symbol, genes)

tbl <- dcast(datp, Hugo_Symbol ~ barcode, value.var="mutation")
tbl <- tbl[order(rev(tbl$Hugo_Symbol)),]

tbl$num <-  apply(tbl[,2:134], 1, function (x) sum(x, na.rm=T))
tbl$cases <-  apply(tbl[,2:134], 1, function (x) length(x[!is.na(x)]))
tbl$freq <-  apply(tbl[,2:134], 1, function (x) (sum(x, na.rm=T)/length(x[!is.na(x)]))*100)

write.table(tbl, "~/brca/brcatcga/analysis/R/mutTable.txt", sep="\t", row.names=F, quote=F)

datp$mutation <- factor(datp$mutation)
#datp$barcode <- factor(datp$barcode, levels=datp$barcode[order(mutSum$imm, mutSum$mutCnt)])

width <- 3125
#ratio <- 461/572
ratio <- 547/600
#pointsize <- 12 x (300/72)
height <- trunc(width*ratio)
#png(paste0("~/brca/brcatcga/analysis/R/top34mutHeatmapv2.png"), width=width, height=height, res=300)
png(paste0("F:\\TNBC TILS\\brcatcga\\analysis\\R\\top34mutHeatmapv2.png"), width=width, height=height, res=300)

# top X TNBC
ggplot(datp, aes(x=barcode, y=Hugo_Symbol, fill=mutation)) + 
  geom_tile() +
  scale_fill_manual(values = c("#FFFFFF","#323D8D"), na.value="#BEBEBE") +
  scale_x_discrete(expand=c(0,1.1)) +
  scale_y_discrete(expand=c(0,0.57)) +
  guides(fill=F) +
  ylab("Genes") +
  xlab("Sample") +
  #ggtitle("Mutations") +
  theme(plot.title = element_text(size=20, hjust = 0.5, face="bold")) +
  theme(text = element_text(size=14, colour="black")) +
  theme(axis.text = element_text(size=12, colour="black")) +
  theme(axis.text.x = element_text(size=8, vjust=0.5, hjust=1, angle=90, colour="black")) +
  theme(axis.text.y = element_text(size=20, colour="black", face="bold")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(panel.grid.major.y=element_blank()) +
  theme(panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.y=element_blank()) +
  theme(panel.grid.minor.x=element_blank()) +
  theme(panel.background=element_rect(fill="#000000")) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())

dev.off()

# sig genes in CD8