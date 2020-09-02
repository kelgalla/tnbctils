



vals <- "abs"
rel <- read.delim("~/brca/brcatcga/analysis/R/cibersort/data/CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

cell_idx <- 5
cibersort <- rel[,c(1,cell_idx)]
cellType <- colnames(cibersort)[2]
colnames(cibersort)[2] <- "sum"

surFile <- "~/brca/brcatcga/analysis/R/tnbc.tils.append.surv.idc.txt"
sur <- read.delim(surFile, header=T, stringsAsFactors=F, na.strings=(""))
#cibersort <- read.delim("~/brca/analysis/R/tnbcGeneExpCiberSort.htseq.txt", header=T, stringsAsFactors=F, na.strings=(""))
sur <- merge(sur, cibersort, by.x="bcr_patient_barcode", by.y="Input.Sample", all.x=T)
#smpl.surv <- Surv(sur$os_days, sur$os_status)~1
sur[sur$pathologic_stage %in% c("Stage I", "Stage IA"),]$pathologic_stage <- "Stage I"
sur[sur$pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB"),]$pathologic_stage <- "Stage II"
sur[sur$pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),]$pathologic_stage <- "Stage III"
#sur <- sur[!(is.na(sur$pathologic_stage)),]
sur$pathologic_stage <- factor(sur$pathologic_stage, levels=c("Stage I", "Stage II", "Stage III", "Stage IV"))

# analysis of cibersort
sur <- sur[!(is.na(sur$sum)),]

sur$cat <- ifelse(sur$sum < 0.1011, "low", ifelse(sur$sum > 0.1726, "high", "medium"))

immQuant <- "low"

# age_at_initial_pathologic_diagnosis
length(sort(sur[sur$cat %in% c("high", "low"),]$age_at_initial_pathologic_diagnosis))

sort(sur[sur$cat==immQuant,]$age_at_initial_pathologic_diagnosis)
length(sort(sur[sur$cat==immQuant,]$age_at_initial_pathologic_diagnosis))

sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(29:39))),]$age_at_initial_pathologic_diagnosis
sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(40:49))),]$age_at_initial_pathologic_diagnosis
sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(50:59))),]$age_at_initial_pathologic_diagnosis
sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(60:69))),]$age_at_initial_pathologic_diagnosis
sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(70:79))),]$age_at_initial_pathologic_diagnosis
sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(80:90))),]$age_at_initial_pathologic_diagnosis

length(sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(29:39))),]$age_at_initial_pathologic_diagnosis)
length(sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(40:49))),]$age_at_initial_pathologic_diagnosis)
length(sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(50:59))),]$age_at_initial_pathologic_diagnosis)
length(sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(60:69))),]$age_at_initial_pathologic_diagnosis)
length(sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(70:79))),]$age_at_initial_pathologic_diagnosis)
length(sur[(sur$cat==immQuant & (sur$age_at_initial_pathologic_diagnosis %in% c(80:90))),]$age_at_initial_pathologic_diagnosis)

# pathologic_stage
length(sort(sur[sur$cat %in% c("high", "low"),]$pathologic_stage, na.last=T))
sum(!is.na(sort(sur[sur$cat %in% c("high", "low"),]$pathologic_stage, na.last=T)))

sort(sur[sur$cat==immQuant,]$pathologic_stage, na.last=T)
length(sort(sur[sur$cat==immQuant,]$pathologic_stage, na.last=T))
sum(!is.na(sort(sur[sur$cat==immQuant,]$pathologic_stage, na.last=T)))

tmp <- sur[(sur$cat==immQuant & sur$pathologic_stage == "Stage I"),]$pathologic_stage
sum(!is.na(tmp))

tmp <- sur[(sur$cat==immQuant & sur$pathologic_stage == "Stage II"),]$pathologic_stage
sum(!is.na(tmp))

tmp <- sur[(sur$cat==immQuant & sur$pathologic_stage == "Stage III"),]$pathologic_stage
sum(!is.na(tmp))

tmp <- sur[(sur$cat==immQuant & sur$pathologic_stage == "Stage IV"),]$pathologic_stage
sum(!is.na(tmp))

tmp <- sur[(sur$cat==immQuant & is.na(sur$pathologic_stage)),]$pathologic_stage
sum(is.na(tmp))

