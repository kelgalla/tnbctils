#---------------------------------------------------------------------
# FILE     : expForCIBERSORT.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-03
# COMMENTS : write out the file to input into CIBERSORT
#---------------------------------------------------------------------

expFile <- 'F:\\TNBC TILS\\brca_metabric_cbioportal\\data_expression_median.txt'
casesFile <- 'F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcCases.txt'
tnbcCasesExpFile <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\tnbcCasesExp.txt"
expForCIBERSORT <- "F:\\TNBC TILS\\tnbctils\\brcametabric\\processedData\\expForCIBERSORT.txt"

cases <- read.delim(casesFile, header=T, stringsAsFactors=F)
exp <- read.delim(expFile, header=T, stringsAsFactors=F)

exp <- subset(exp, select = -c(Entrez_Gene_Id)) # remove Entrez column

exp <- exp[,colnames(exp) %in% c("Hugo_Symbol",cases$PATIENT_ID)] # keep only TNBC cases, 196 now

# reverse log space for cibersort
exp[,2:length(exp)] <- 2^exp[,2:length(exp)]

# save cases with expression
casesExp <- colnames(exp[,-1])

write.table(data.frame(PATIENT_ID=casesExp), tnbcCasesExpFile, sep="\t", quote=F, row.names=F)

# save CIBERSORT input file
write.table(exp, expForCIBERSORT, sep="\t", quote=F, row.names=F)
