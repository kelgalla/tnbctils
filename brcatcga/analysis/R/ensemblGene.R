#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library(biomaRt)

htseq <- read.delim("~/brca/analysis/bin/tnbcGeneExp.fpkm.txt", header=T, stringsAsFactors=F)
names(htseq)[names(htseq) == 'ensembl_gene_id'] <- "ensembl_version"

# remove __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique
#toRm <- htseq[(htseq$ensembl_version %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")),]
#htseq <- htseq[!(htseq$ensembl_version %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")),]

# remove version
htseq$ensembl_gene_id <- gsub("\\.[0-9]+", "", htseq$ensembl_version)

# move new col to front
col_idx <- grep("^ensembl_gene_id$", names(htseq), perl=T)
htseq <- htseq[, c(col_idx, (1:ncol(htseq))[-col_idx])]

#ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
hgnc <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','hgnc_id', 'external_gene_name'),filters = 'ensembl_gene_id', values = htseq$ensembl_gene_id, mart = ensembl)
hgnc[hgnc$hgnc_symbol=="",]$hgnc_symbol <- NA
hgnc[hgnc$hgnc_id=="",]$hgnc_id <- NA

# there are some duplicated genes
hgnc$ensembl_gene_id[duplicated(hgnc$ensembl_gene_id) | duplicated(hgnc$ensembl_gene_id, fromLast=TRUE)]

temp2 <- hgnc[!is.na(hgnc$hgnc_symbol),]
temp2$ensembl_gene_id[duplicated(temp2$ensembl_gene_id) | duplicated(temp2$ensembl_gene_id, fromLast=TRUE)]
temp2$hgnc_symbol[duplicated(temp2$hgnc_symbol) | duplicated(temp2$hgnc_symbol, fromLast=TRUE)]
temp3 <- temp2[temp2$ensembl_gene_id=="ENSG00000278233",]
temp3 <- temp2[temp2$hgnc_symbol=="SNORD38B",]

htseqFin <- merge(htseq, hgnc, by="ensembl_gene_id", all.x=T)

# move new cols to front
egi_idx <- grep("^ensembl_gene_id$", names(htseqFin), perl=T)
egiv_idx <- grep("^ensembl_version$", names(htseqFin), perl=T)
hgncSym_idx <- grep("^hgnc_symbol$", names(htseqFin), perl=T)
hgncID_idx <- grep("^hgnc_id$", names(htseqFin), perl=T)
geneName_idx <- grep("^external_gene_name$", names(htseqFin), perl=T)
htseqFin <- htseqFin[, c(egi_idx, egiv_idx, hgncSym_idx, hgncID_idx, geneName_idx, (1:ncol(htseq))[-c(egi_idx, egiv_idx, hgncSym_idx, hgncID_idx, geneName_idx)])]

# in the output file, it will be hgnc_symbol and samples only, so remove na rows, it is ok that some symbols are still redundant
htseqOut <- htseqFin[!(is.na(htseqFin$hgnc_symbol)),]
hgncSym_idx <- grep("^hgnc_symbol$", names(htseqOut), perl=T)
htseqOut <- htseqOut[,c(hgncSym_idx, 6:dim(htseqOut)[2])]
write.table(htseqOut, "~/brca/analysis/R/tnbcGeneExpSymbol.fpkm.txt", sep="\t", quote=F, row.names=F)


temp <- htseqFin[!(is.na(htseqFin$hgnc_symbol)),]
temp$hgnc_symbol[duplicated(temp$hgnc_symbol) | duplicated(temp$hgnc_symbol, fromLast=TRUE)]
temp$hgnc_id[duplicated(temp$hgnc_id) | duplicated(temp$hgnc_id, fromLast=TRUE)]

temp <- htseqFin[!(is.na(htseqFin$hgnc_symbol)),]
temp <- temp[temp$hgnc_symbol == "EMG1",]

temp <- temp[is.na(temp$hgnc_symbol),]
temp <- temp[!(temp$hgnc_symbol == temp$external_gene_name),]

temp <- htseqFin[is.na(htseqFin$external_gene_name),]
sum(is.na(htseqFin$external_gene_name))
sum(is.na(htseqFin$hgnc_symbol))

att <- listAttributes(ensembl)

