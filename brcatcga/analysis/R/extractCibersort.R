


rel <- read.delim("~/brca/analysis/ANALYSIS/CIBERSORT.Output_Job4.txt", header=T, stringsAsFactors=F)
rel$sum <- rel$T.cells.CD4.memory.activated
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)
out <- rel[,c(1,27)]
write.table(out, "~/brca/analysis/R/tnbcGeneExpCiberSort.htseq.txt", sep="\t", quote=F, row.names=F)

abs <- read.delim("~/brca/analysis/ANALYSIS/CIBERSORT.Output_Abs_Job4.txt", header=T, stringsAsFactors=F)
#cells <- abs[,c(2:23)]
#abs$sum <- rowSums(cells)
abs$sum <- abs$T.cells.CD4.memory.resting

#mean(abs$sum)
#median(abs$sum)
#grp1 <- abs[abs$sum < 26.7,]
#grp2 <- abs[abs$sum > 26.7,]

abs$Input.Sample <- gsub("\\.", "-", abs$Input.Sample)
out <- abs[,c(1,27)]
write.table(out, "~/brca/analysis/R/tnbcGeneExpCiberSort.htseq.txt", sep="\t", quote=F, row.names=F)
