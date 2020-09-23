#---------------------------------------------------------------------
# FILE     : ttMutAnalysis.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-22
# COMMENTS : Explore if mutations tend to occur in the same subtype
#---------------------------------------------------------------------

dir <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R"

#------- CD8 T cells
mutSigDiff <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD8.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)

#------- CD4 memory activated T cells
mutSigDiff <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.abs.T.cells.CD4.memory.activated.0.25sd.mutSigDiff.txt", header=T, stringsAsFactors=F)

#------- CD8 CD4
mutSigDiff <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\os.clusterGrp.mutSigDiff.txt", header=T, stringsAsFactors=F)

genes <- mutSigDiff$Hugo_Symbol

# import immune status of sample for sorting
#------- CD8 T cells
pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD8.0.25sd.txt"), header=T, stringsAsFactors=F)

#------- CD4 memory activated T cells
pts <- read.delim(paste0(dir, "\\os.abs.T.cells.CD4.memory.activated.0.25sd.txt"), header=T, stringsAsFactors=F)

#------- CD8 CD4
pts <- read.delim(paste0(dir, "\\clusterGrp.txt"), header=T, stringsAsFactors=F)
pts <- pts[pts$clusterGrp %in% c("bothPos", "bothNeg"),]
pts[pts$clusterGrp=="bothPos",]$clusterGrp <- "high"
pts[pts$clusterGrp=="bothNeg",]$clusterGrp <- "low"
colnames(pts)[colnames(pts)=="clusterGrp"] <- "imm"

pts$barcode <- gsub("-", "\\.", pts$barcode)

# import mutation count of sample for sorting
mutCnt <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\mutCnt.txt", header=T, stringsAsFactors=F)

# combine immune and mutation data
mutSum <- merge(mutCnt, pts, by.x="patient", by.y="barcode")
mutSum <- mutSum[!is.na(mutSum$imm),] # remove imm NAs
mutSum <- mutSum[!is.na(mutSum$mutCnt),] # remove mutCnt NAs
mutSum <- mutSum[order(mutSum$imm, mutSum$mutCnt),]

# get non-summarized mutation data for graphing
datNew <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\mut.txt", header=T, stringsAsFactors=F)
datNew$barcode <- gsub("-", "\\.", datNew$barcode)
datp <- datNew

datp <- datp[datp$Hugo_Symbol %in% genes,] # get subset of genes to graph
datp$Hugo_Symbol <- factor(datp$Hugo_Symbol, rev(genes))

datp <- datp[datp$barcode %in% mutSum$patient,] # get subset of samples to graph
datp$barcode <- factor(datp$barcode, levels=mutSum$patient) # put in desired order

# high low grps
datp$group <- ifelse(datp$barcode %in% mutSum[mutSum$imm=="high",]$patient, "high" , "low")

ttFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbctype\\sexprs\\bac566ef-2fce-40d5-a6c0-d8678c8a0cf3_result.csv"

tt <- read.csv(ttFile, header=T, stringsAsFactors=F, na.strings=(""))
colnames(tt)[1] <- "barcode"
tt <- tt[,c(1,2)]

datp <- merge(datp, tt, by="barcode", all.x=T)

tmp <- datp[datp$mutation=="1",]
tmp <- tmp[rev(order(tmp$Hugo_Symbol)),]

#------- CD8 T cells
write.table(tmp, paste0(dir, "\\os.abs.T.cells.CD8.0.25sd.mutSubtype.txt"), row.names=F, sep="\t", quote=F)

#------- CD8 CD4
write.table(tmp, paste0(dir, "\\clusterGrp.mutSubtype.txt"), row.names=F, sep="\t", quote=F)
