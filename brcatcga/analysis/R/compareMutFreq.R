#---------------------------------------------------------------------
# FILE     : compareMutFreq.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2018-02-26
# COMMENTS : this file will take a mutation table with high and low 
#            counts, and determine if there is a statistical 
#            difference between the counts for each gene
#---------------------------------------------------------------------

#---- CD8 T Cells
tbl <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.mutTable.txt", header=T, stringsAsFactors=F)
#tbl <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.mediansd.mutTable.txt", header=T, stringsAsFactors=F)

#---- CD4 memory activated T cells
#tbl <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutTable.txt", header=T, stringsAsFactors=F)
tbl <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.mediansd.mutTable.txt", header=T, stringsAsFactors=F)

#---- CD8 CD4
tbl <- read.delim("F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.clusterGrp.mutTable.txt", header=T, stringsAsFactors=F)

end <- tbl[1:8453,c(1,135:143)]
bottom <- tail(end)

hcases <- tbl$highCases[1]
lcases <- tbl$lowCases[1]

genes <- tbl[1:8453,1]

tbl$pval <- NA
#tbl[1:5,135:144]

for (i in genes) {
  
  #i <- "TP53"
  hmut <- ifelse(tbl[tbl$Hugo_Symbol==i,]$highNum==0, 0, tbl[tbl$Hugo_Symbol==i,]$highNum)
  lmut <- ifelse(tbl[tbl$Hugo_Symbol==i,]$lowNum==0, 0, tbl[tbl$Hugo_Symbol==i,]$lowNum)
  
  hnomut <- hcases - hmut
  lnomut <- lcases - lmut
  
  square <- matrix(c(hmut, lmut, hnomut, lnomut), nrow=2)
  
  pval <- fisher.test(square)$p.value
  
  tbl[tbl$Hugo_Symbol==i,]$pval <- pval
  
  if (pval < 0.05){
    print(paste0(i, "\t", pval)) 
  }
}

end <- tbl[1:8453,c(1,135:144)]

padj <- p.adjust(end$pval, method="BH")
end$padj <- padj

goi <- end[end$padj < 0.05,]

goi <- end[end$pval < 0.05,]
goi <- goi[order(goi$pval),]

#------- CD8 T cells
#write.table(goi, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.mutSigDiff.txt", sep="\t", row.names=F, quote=F)
#write.table(goi, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.mediansd.mutSigDiff.txt", sep="\t", row.names=F, quote=F)

#------- CD4 memory activated T cells
#write.table(goi, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", sep="\t", row.names=F, quote=F)
write.table(goi, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.mediansd.mutSigDiff.txt", sep="\t", row.names=F, quote=F)

#------- CD8 CD4
#write.table(goi, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", sep="\t", row.names=F, quote=F)
write.table(goi, "F:\\TNBC TILS\\brcatcga\\analysis\\R\\os.clusterGrp.mutSigDiff.txt", sep="\t", row.names=F, quote=F)

mutSigDiffcd8 <- read.delim("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD8.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
#mutSigDiffcd8 <- read.delim("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD8.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
mutSigDiffcd4 <- read.delim("~/brca/brcatcga/analysis/R/os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)
