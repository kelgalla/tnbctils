#---------------------------------------------------------------------
# FILE     : clinicalTNBC.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-03
# COMMENTS : identify TNBC cases with associated survival data, 
#            write it out to files
#---------------------------------------------------------------------

#-------
# if using synapse version of data, need to convert illumina probes to hugo gene wheras cbioportal has already done this
#library("illuminaHumanv4.db")
#synapse <- 'F:\\TNBC TILS\\brca_metabric_synapse\\EGAD00010000210\\discovery_ExpressionMatrix.txt'
#synapse <- 'F:\\TNBC TILS\\brca_metabric_synapse\\EGAD00010000211\\validation_ExpressionMatrix.txt'
#syndata <- read.delim(synapse, header=T, stringsAsFactors=F, sep=' ')
#probeID <- c("ILMN_1802380")
#data.frame(Gene=unlist(mget(x = probeID,envir = illuminaHumanv4SYMBOL)))
#--------

pt <- 'F:\\TNBC TILS\\brca_metabric_cbioportal\\data_clinical_patient.txt'
smpl <- 'F:\\TNBC TILS\\brca_metabric_cbioportal\\data_clinical_sample.txt'

tnbcSurv <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcSurv.txt"
tnbcCases <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcCases.txt"

# filter TNBC columns in patient file
ptdata <- read.delim(pt, header=T, stringsAsFactors=F, comment.char="#")
pttnbc <- ptdata[!is.na(ptdata$ER_IHC) & ptdata$ER_IHC=="Negative",]
pttnbc <- pttnbc[!is.na(pttnbc$HER2_SNP6) & !pttnbc$HER2_SNP6=="GAIN",] # though some called negative in smpl file, don't want gain by snp
pttnbc <- pttnbc[!is.na(pttnbc$HISTOLOGICAL_SUBTYPE) & pttnbc$HISTOLOGICAL_SUBTYPE=="Ductal/NST",]

# filter to TNBC columns in sample file
smpldata <- read.delim(smpl, header=T, stringsAsFactors=F, comment.char="#")
smpltnbc <- smpldata[!is.na(smpldata$ER_STATUS) & smpldata$ER_STATUS=="Negative",]
smpltnbc <- smpltnbc[!is.na(smpltnbc$PR_STATUS) & smpltnbc$PR_STATUS=="Negative",]
smpltnbc <- smpltnbc[!is.na(smpltnbc$HER2_STATUS) & smpltnbc$HER2_STATUS=="Negative",]
smpltnbc <- smpltnbc[!is.na(smpltnbc$CANCER_TYPE_DETAILED) & smpltnbc$CANCER_TYPE_DETAILED=="Breast Invasive Ductal Carcinoma",]

# no case of multiple samples per patient
#all(smpltnbc$PATIENT_ID==smpltnbc$SAMPLE_ID)

# merge info from sample and patient files
tnbc <- merge(pttnbc, smpltnbc) # 199 cases

# write out survival data and case id list
write.table(tnbc, tnbcSurv, sep="\t", quote=F, row.names=F)

cases <- sapply(data.frame(PATIENT_ID=tnbc$PATIENT_ID), function(x) gsub("-", ".", x))

write.table(cases, tnbcCases, sep="\t", quote=F, row.names=F)






