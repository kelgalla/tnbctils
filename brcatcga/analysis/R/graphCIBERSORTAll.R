#---------------------------------------------------------------------
# FILE     : graphCIBERSORTAll.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-16
# COMMENTS : This script takes in the cibersort results, and graphs 
#            a clustering heatmap for all immune types
#---------------------------------------------------------------------

library(gplots)
library(Heatplus)
library(dendextend)

#-----------------------
# All immune cells
cibersort <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)

dirSave <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\"

cm <- as.matrix(cibersort[,c(2:5,7:23)]) # column 6 is all 0's
rownames(cm) <- cibersort[,1]

exprs <- cm
sexprs <- scale(exprs)

my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)



col_breaks <- c(seq(-1, -0.25, length=100),
                seq(-0.245, 0.245, length=100),
                seq(0.25, 4.4, length=100))

hold <- heatmap.2(sexprs, Rowv=T, Colv=NA, dendrogram="row", lwid=c(1,2),
                  na.color="#D3D3D3", col=my_palette, trace="none",
                  scale="none", cexRow=1, cexCol=1, labRow=NA, labCol=NA,
                  density.info="none", symbreaks=T, symkey=F, key=T, margins=c(4,4), breaks=col_breaks)

#------------
# cluster analysis

clusterGrp <- data.frame(barcode=rownames(sexprs), clusterGrp="tbd", stringsAsFactors=F)

# first cut
hr <- as.dendrogram(hclust(dist(sexprs)))
Rowv <- rowMeans(sexprs)
hr <- reorder(hr, Rowv)
hrclust <- as.hclust(hr)

mycl <- cutree(hrclust, k=11) # 11 clusters
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); 
mycolhc <- mycolhc[as.vector(mycl)]
names(mycolhc) <- names(mycl)

width <- 1750
#ratio <- 547/600
ratio <- 0.9 # height/width
height <- trunc(width*ratio)
png(paste0(dirSave, "CIBERSORT.Output_Abs_fpkm.heatmapbycolorscluster.png"), 
    width=width, height=height, res=300)

hm <- heatmap.2(sexprs, Rowv=hr, Colv=NULL, dendrogram="row", lwid=c(1,2),
                na.color="#D3D3D3", col=my_palette, trace="none",
                scale="none", cexRow=1, cexCol=1, labRow=NA, srtCol = 45,
                density.info="none", symbreaks=T, symkey=F, key=T, margins=c(10,4), breaks=col_breaks, RowSideColors=mycolhc)

dev.off()

x11(height=10, width=2); 
barplot(rep(10, max(mycl)), col=unique(mycolhc[hrclust$labels[hrclust$order]]), horiz=T, names=unique(mycl[hrclust$order]))
