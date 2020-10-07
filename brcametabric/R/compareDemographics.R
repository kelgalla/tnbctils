#---------------------------------------------------------------------
# FILE     : compareDemographics.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-10-05
# COMMENTS : file for comparing tcga and metabric demographics
#---------------------------------------------------------------------


# tcga
tcgasurFile <- "F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\tnbc.tils.append2.surv.idc.txt"
tcgasur <- read.delim(tcgasurFile, header=T, stringsAsFactors=F, na.strings=(""))

# metabric
metabricsurFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcSurv.txt"
metabricsur <- read.delim(metabricsurFile, header=T, stringsAsFactors=F, na.strings=(""))
metabricsur$OS_STATUS2 <- ifelse(metabricsur$OS_STATUS == "0:LIVING", 0, 1)

# age
wilcox.test(tcgasur$age_at_initial_pathologic_diagnosis, metabricsur$AGE_AT_DIAGNOSIS)

# stage
t <- matrix(c(110, 20, 3, 129, 12, 58), 3, 2, dimnames=list(c("stage I-II", "stage III-IV", "na"), c("tcga", "metabric")))
fisher.test(t)

t <- matrix(c(110, 20, 129, 12), 2, 2, dimnames=list(c("stage I-II", "stage III-IV"), c("tcga", "metabric")))
fisher.test(t)

# menopause status
t <- matrix(c(35,82,16,82,117,0), 3, 2, dimnames=list(c("pre", "post", "na"), c("tcga", "metabric")))
fisher.test(t)

t <- matrix(c(35,82,82,117), 2, 2, dimnames=list(c("pre", "post"), c("tcga", "metabric")))
fisher.test(t)

# chemotherapy
tcgachemo <- !is.na(tcgasur$drug_name)
tcgachemoY <- sum(tcgachemo, na.rm=T)
tcgachemoN <- length(tcgachemo) - tcgachemoY

metabricchemoY <- nrow(metabricsur[metabricsur$CHEMOTHERAPY=="YES",])
metabricchemoN <- nrow(metabricsur[metabricsur$CHEMOTHERAPY=="NO",])

t <- matrix(c(103,30,117,82), 2, 2, dimnames=list(c("chemo", "no chemo"), c("tcga", "metabric")))
fisher.test(t)

# LN status
tcgaLNNeg <- nrow(tcgasur[tcgasur$pathologic_N=="N0" | tcgasur$pathologic_N=="N0 (i-)" | tcgasur$pathologic_N=="N0 (i+)",])
tcgaLNPos <- nrow(tcgasur) - tcgaLNNeg

metabricsur$LYMPH_NODES_EXAMINED_POSITIVE <- as.numeric(metabricsur$LYMPH_NODES_EXAMINED_POSITIVE)
metabricLNNeg <- sum(metabricsur$LYMPH_NODES_EXAMINED_POSITIVE == 0, na.rm=T)
metabricLNPos <- sum(metabricsur$LYMPH_NODES_EXAMINED_POSITIVE > 0 , na.rm=T)

t <- matrix(c(49,84,0,100,96,3), 3, 2, dimnames=list(c("LNPos", "LNNeg", "NA"), c("tcga", "metabric")))
fisher.test(t)

t <- matrix(c(49,84,100,96), 2, 2, dimnames=list(c("LNPos", "LNNeg"), c("tcga", "metabric")))
fisher.test(t)
