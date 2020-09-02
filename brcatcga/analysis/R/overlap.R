

dir <- "~/brca/dataLegacyReorg/"
badveFile <- paste0(dir, "TNBCfromBadve.txt")
meFile <- paste0(dir,"tnbc.txt")

badve <- read.delim(badveFile, header=T, stringsAsFactors=F)
badve <- data.frame("bcr_patient_barcode" = badve$bcr_patient_barcode, "badve" = 1)
me <- read.delim(meFile, header=T, stringsAsFactors=F)
me <- data.frame("bcr_patient_barcode" = me$bcr_patient_barcode, "me" = 1)

match <- merge(badve, me, by.x="bcr_patient_barcode", all = T)
