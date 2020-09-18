#---------------------------------------------------------------------
# FILE     : graphCIBERSORT.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2018-02-27
# COMMENTS : This script takes in the cibersort results, and graphs 
#            heatmaps of the columns of interest
#---------------------------------------------------------------------

library(gplots)
library(Heatplus)
library(dendextend)

#-----------------------
# T Cells CD8
cibersort <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)

cibersort <- cibersort[,c("Input.Sample", "T.cells.CD8")]
#cibersort <- cibersort[order(cibersort$T.cells.CD8),]

cm <- as.matrix(cibersort[,c(2,2)]) 
rownames(cm)<- cibersort[,1]

exprs <- cm

#my_palette <- colorRampPalette(c("blue", "grey", "red"))(n = 299)
my_palette <- colorRampPalette(c("#F8766D", "black", "#00BFC4"))(n = 299)

col_breaks = c(seq(0, 0.101, length=100),
               seq(0.1011,0.1725, length=100),
               seq(0.1726, 0.6908, length=100))

#my_palette <- colorRampPalette(c("blue", "red"))(n = 300)

png("~/brca/brcatcga/analysis/R/CIBERSORT.Output_Abs_fpkm.T.cells.CD8.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300)

#lmat = rbind(c(0,4),c(2,1),c(0,3))
#lwid = c(1.5,4)
#lhei = c(1.5,4,1.5)
#margins=c(10,11)
#srtCol=45

hold <- heatmap.2(exprs, Rowv=NA, Colv=NA, dendrogram="none", lwid=c(1,2),
                  na.color="#D3D3D3", col=my_palette, trace="none", 
                  scale="none", cexRow=1, cexCol=1, labCol=NA, labRow=NA,
                  density.info="none", breaks=col_breaks, symbreaks=F, key=T, margins=c(2,7))

#                  tracecol="#303030", sepwidth=c(0.4,0.4), sepcolor="white", rowsep=1:nrow(exprs))

dev.off()


#-----------------------
# T Cells CD4 Memory Activated
cibersort2 <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)

cibersort2 <- cibersort2[,c("Input.Sample", "T.cells.CD4.memory.activated")]
cibersort2 <- cibersort2[order(cibersort2$T.cells.CD4.memory.activated),]

cm <- as.matrix(cibersort2[,c(2,2)]) 
rownames(cm)<- cibersort2[,1]

exprs <- cm

#my_palette <- colorRampPalette(c("blue", "grey", "red"))(n = 299)
my_palette <- colorRampPalette(c("#F8766D", "black", "#00BFC4"))(n = 299)

col_breaks = c(seq(0, 0.0347, length=100),
               seq(0.0348,0.0645, length=100),
               seq(0.0646, 0.3079, length=100))

#my_palette <- colorRampPalette(c("blue", "red"))(n = 300)

png("~/brca/brcatcga/analysis/R/CIBERSORT.Output_Abs_fpkm.T.cells.CD4.memory.activated.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300)

#lmat = rbind(c(0,4),c(2,1),c(0,3))
#lwid = c(1.5,4)
#lhei = c(1.5,4,1.5)
#margins=c(10,11)
#srtCol=45

hold <- heatmap.2(exprs, Rowv=NA, Colv=NA, dendrogram="none", lwid=c(1,2),
                  na.color="#D3D3D3", col=my_palette, trace="none", 
                  scale="none", cexRow=1, cexCol=1, labCol=NA, labRow=NA,
                  density.info="none", breaks=col_breaks, symbreaks=F, key=T, margins=c(2,7))

#                  tracecol="#303030", sepwidth=c(0.4,0.4), sepcolor="white", rowsep=1:nrow(exprs))

dev.off()

#-----------------------
# T Cells CD8 and T Cells CD4
cibersort <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)

cibersort <- cibersort[,c("Input.Sample", "T.cells.CD8", "T.cells.CD4.memory.activated")]
cibersort <- cibersort[order(cibersort$T.cells.CD8),]

cm <- as.matrix(cibersort[,c(2,3)])
rownames(cm) <- cibersort[,1]

exprs <- cm
sexprs <- scale(exprs)

#my_palette <- colorRampPalette(c("#F8766D", "black", "#00BFC4"))(n = 299)
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)

#col_breaks <- c(seq(min(c(sexprs[,1],c(sexprs[,2]))), -0.25, length=150),
#                seq(0.25,max(c(sexprs[,1], sexprs[,2])), length=150))

col_breaks <- c(seq(-1, -0.25, length=100),
                seq(-0.245, 0.245, length=100),
                seq(0.25, 4.4, length=100))


#col_breaks = c(seq(0, 0.101, length=100),
#               seq(0.1011,0.1725, length=100),
#               seq(0.1726, 0.6908, length=100))

#my_palette <- colorRampPalette(c("blue", "red"))(n = 300)

#png("~/brca/brcatcga/analysis/R/CIBERSORT.Output_Abs_fpkm.T.cells.CD8.png",    # create PNG for the heat map        
#    width = 5*300,        # 5 x 300 pixels
#    height = 5*300,
#    res = 300)

#lmat = rbind(c(0,4),c(2,1),c(0,3))
#lwid = c(1.5,4)
#lhei = c(1.5,4,1.5)
#margins=c(10,11)
#srtCol=45
#col=g2r.colors(256, min.tinge=0.33)

#pear = function(x) as.dist(1-cor(t(x), method="pearson"))
#avg = function(x) hclust(x, method="average")
#corrdist = function(x) dist(x, method="euclidian")
#hclust.avl = function(x) hclust(x, method="complete")

#hr <- as.dendrogram(avg(pear(exprs)))
#Rowv <- rowMeans(sexprs)
#hr <- reorder(hr, Rowv, max) # default is sum

#hr <- as.dendrogram(hclust.avl(corrdist(sexprs)))
#Colv <- rowMeans(t(sexprs))
#hc <- rev(reorder(hc, Colv)) # default is sum

width <- 1750
#ratio <- 547/600
ratio <- 0.9 # height/width
height <- trunc(width*ratio)
#png("~/brca/brcatcga/analysis/R/CIBERSORT.Output_Abs_fpkm.T.cells.CD8.T.cells.CD4.memory.activated.png", width=width, height=height, res=300)
png("F:\\TNBC TILS\\brcatcga\\analysis\\R\\CIBERSORT.Output_Abs_fpkm.T.cells.CD8.T.cells.CD4.memory.activated.heatmapbycolors.png", 
    width=width, height=height, res=300)

#hold <- heatmap.2(sexprs, Rowv=T, Colv=NA, dendrogram="row", lwid=c(1,2),
#                  na.color="#D3D3D3", col=g2r.colors(256, min.tinge=0.33), trace="none", symkey=T,
#                  scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
#                  density.info="none", symbreaks=T, key=T, margins=c(4,4))

#par(cex=1.5)
hold <- heatmap.2(sexprs, Rowv=T, Colv=NA, dendrogram="row", lwid=c(1,2),
                  na.color="#D3D3D3", col=my_palette, trace="none",
                  scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
                  density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks)


# col_breaks = c(seq(0, 0.101, length=150),
#                seq(0.1011,0.6908, length=150))
# 
# col_breaks = c(seq(0, 0.101, length=100),
#                seq(0.1011,0.1725, length=100),
#                seq(0.1726, 0.6908, length=100))
# 
# hold <- heatmap.2(exprs, Rowv=T, Colv=NA, dendrogram="row", lwid=c(1,2),
#                   na.color="#D3D3D3", col=g2r.colors(299, min.tinge=0.33), trace="none", symkey=F,
#                   scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
#                   density.info="none", symbreaks=F, key=T, margins=c(4,4), breaks=col_breaks)
# 
# #                  tracecol="#303030", sepwidth=c(0.4,0.4), sepcolor="white", rowsep=1:nrow(exprs))

dev.off()

#------------
# cluster analysis

clusterGrp <- data.frame(barcode=rownames(sexprs), clusterGrp="tbd", stringsAsFactors=F)

# first cut
hr <- as.dendrogram(hclust(dist(sexprs)))
Rowv <- rowMeans(sexprs)
hr <- reorder(hr, Rowv)
hrclust <- as.hclust(hr)

mycl <- cutree(hrclust, k=11)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); 
mycolhc <- mycolhc[as.vector(mycl)]
names(mycolhc) <- names(mycl)

# 1 #FF9900FF
# 2 #EBFF00FF
# 3 #70FF00FF
# 4 #00FF0AFF
# 5 #00FF85FF
# 6 #00FFFFFF
# 7 #0085FFFF
# 8 #000AFFFF
# 9 #7000FFFF
# 10 #EB00FFFF
# 11 #FF0099FF

hm <- heatmap.2(sexprs, Rowv=hr, Colv=NULL, dendrogram="row", lwid=c(1,2),
          na.color="#D3D3D3", col=my_palette, trace="none",
          scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
          density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks, RowSideColors=mycolhc)

x11(height=10, width=2); 
barplot(rep(10, max(mycl)), col=unique(mycolhc[hrclust$labels[hrclust$order]]), horiz=T, names=unique(mycl[hrclust$order]))

# grouping colors
mycl2 <- mycl
mycolhc2 <- mycolhc

mycolhc2[mycl2%in%c(8, 11, 5, 4, 6, 10, 9)] <- rep("#EB00FFFF",length(mycolhc2[mycl2%in%c(8, 11, 5, 4, 6, 10, 9)]))
mycolhc2[mycl2%in%c(7, 3)] <- rep("#00FF0AFF",length(mycolhc2[mycl2%in%c(7,3)]))
mycolhc2[mycl2%in%c(2)] <- rep("#FF0099FF",length(mycolhc2[mycl2%in%c(2)]))

width <- 1800
#ratio <- 547/600
ratio <- 0.9 # height/width
height <- trunc(width*ratio)
png("F:\\TNBC TILS\\brcatcga\\tnbctils\\analysis\\R\\CIBERSORT.Output_Abs_fpkm.T.cells.CD8.T.cells.CD4.memory.activated.heatmapbycolorscluster.png", width=width, height=height, res=300)

hm <- heatmap.2(sexprs, Rowv=hr, Colv=NULL, dendrogram="row", lwid=c(1,2),
                na.color="#D3D3D3", col=my_palette, trace="none",
                scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
                density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks, RowSideColors=mycolhc2)

dev.off()

# both neg, cluster 1
clusterKeep <- c(1)
clusterPrune <- c(3,2,7,9,10,6,4,5,11,8)
toKeep <- names(mycl[mycl%in%clusterKeep])
toPrune <- names(mycl[mycl%in%clusterPrune])
sub <- sexprs[toKeep,]
hrsub <- prune(hr, toPrune)

heatmap.2(sub, Rowv=hrsub, Colv=NA, dendrogram="row", lwid=c(1,2),
          na.color="#D3D3D3", col=my_palette, trace="none",
          scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
          density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks, RowSideColors=mycolhc[mycl%in%clusterKeep])

clusterGrp[(clusterGrp$barcode %in% toKeep),]$clusterGrp <- "bothNeg"
nrow(clusterGrp[clusterGrp$clusterGrp=="bothNeg",])

# figure out what to do with the rest
clusterKeep <- c(3,2,7,9,10,6,4,5,11,8);
clusterPrune <- c(1)
toKeep <- names(mycl[mycl%in%clusterKeep])
toPrune <- names(mycl[mycl%in%clusterPrune])
sub <- sexprs[toKeep,]
hrsub <- prune(hr, toPrune)

heatmap.2(sub, Rowv=hrsub, Colv=NA, dendrogram="row", lwid=c(1,2),
          na.color="#D3D3D3", col=my_palette, trace="none",
          scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
          density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks, RowSideColors=mycolhc[mycl%in%clusterKeep])

# both up, c(8, 11, 5, 4, 6, 10, 9)
clusterKeep <- c(8, 11, 5, 4, 6, 10, 9);
clusterPrune <- c(7, 2, 3, 1)
toKeep <- names(mycl[mycl%in%clusterKeep])
toPrune <- names(mycl[mycl%in%clusterPrune])
sub <- sexprs[toKeep,]
hrsub <- prune(hr, toPrune)

heatmap.2(sub, Rowv=hrsub, Colv=NA, dendrogram="row", lwid=c(1,2),
          na.color="#D3D3D3", col=my_palette, trace="none",
          scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
          density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks, RowSideColors=mycolhc[mycl%in%clusterKeep])

clusterGrp[(clusterGrp$barcode %in% toKeep),]$clusterGrp <- "bothPos"
nrow(clusterGrp[clusterGrp$clusterGrp=="bothPos",])

# cd4 up, cd8 down c(2)
clusterKeep <- c(2);
clusterPrune <- c(8, 11, 5, 4, 6, 10, 9, 7, 3, 1)
toKeep <- names(mycl[mycl%in%clusterKeep])
toPrune <- names(mycl[mycl%in%clusterPrune])
sub <- sexprs[toKeep,]
hrsub <- prune(hr, toPrune)

heatmap.2(sub, Rowv=hrsub, Colv=NA, dendrogram="row", lwid=c(1,2),
          na.color="#D3D3D3", col=my_palette, trace="none",
          scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
          density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks, RowSideColors=mycolhc[mycl%in%clusterKeep])

clusterGrp[(clusterGrp$barcode %in% toKeep),]$clusterGrp <- "cd4upcd8down"
nrow(clusterGrp[clusterGrp$clusterGrp=="cd4upcd8down",])

# cd8up, cd4 down c(7, 3)
clusterKeep <- c(7, 3);
clusterPrune <- c(8, 11, 5, 4, 6, 10, 9, 2, 1)
toKeep <- names(mycl[mycl%in%clusterKeep])
toPrune <- names(mycl[mycl%in%clusterPrune])
sub <- sexprs[toKeep,]
hrsub <- prune(hr, toPrune)

heatmap.2(sub, Rowv=hrsub, Colv=NA, dendrogram="row", lwid=c(1,2),
          na.color="#D3D3D3", col=my_palette, trace="none",
          scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
          density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks, RowSideColors=mycolhc[mycl%in%clusterKeep])

clusterGrp[(clusterGrp$barcode %in% toKeep),]$clusterGrp <- "cd8upcd4down"
nrow(clusterGrp[clusterGrp$clusterGrp=="cd8upcd4down",])

write.table(clusterGrp, "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\clusterGrp.txt", sep="\t", row.names=F, quote=F)
