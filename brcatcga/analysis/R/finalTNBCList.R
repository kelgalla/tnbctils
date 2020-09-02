#---------------------------------------------------------------------
# FILE     : finalTNBCList.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-08-24
# COMMENTS : read in the R/tnbc.tils.append.surv.txt and filter it to
#            a final list to go into cibersort
#---------------------------------------------------------------------

tnbc <- read.delim("~/brca/analysis/R/tnbc.tils.append.surv.txt", header=T, stringsAsFactors=F)

keep <- tnbc[tnbc$histological_type == "Infiltrating Ductal Carcinoma",]

dontkeep <- tnbc[tnbc$histological_type != "Infiltrating Ductal Carcinoma",]

write.table(keep, "~/brca/analysis/R/tnbc.tils.append.surv.idc.txt", sep="\t", quote=F, row.names=F)
