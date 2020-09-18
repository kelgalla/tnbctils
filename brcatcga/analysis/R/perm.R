#---------------------------------------------------------------------
# FILE     : perm.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2020-09-18
# COMMENTS : This file will check to see if any of the cibersort 
#            returned values are effected by using 1000 permutations 
#            vs. 100 permutations
#---------------------------------------------------------------------


#rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

#rel1000 <- read.delim("F:\\TNBC TILS\\tnbctils\\brcatcga\\analysis\\R\\cibersort\\data1000\\CIBERSORT.Output_Abs_fpkm.txt", header=T, stringsAsFactors=F)
rel1000 <- read.delim("F:\\TNBC TILS\\tnbctils\\brcametabric\\cibersort\\data1000\\CIBERSORT.Output_Abs_intensity.txt", header=T, stringsAsFactors=F)
rel1000$Input.Sample <- gsub("\\.", "-", rel$Input.Sample)

print(rel[1,])
print(rel1000[1,])

for (row in 1:nrow(rel)){
  if (!identical(rel[row,c(1:23,25,26)], rel1000[row,c(1:23,25,26)])){
    print(rel[row,])
    print(rel1000[row,])
  }
}

# pvalue row
for (row in 1:nrow(rel)){
  if (!identical(rel[row,24], rel1000[row,24])){
      print(rel[row,c(1,24)])
      print(rel1000[row,c(1,24)])
  }
}

relp <- rel[,c(1,24)]
relp1000 <- rel1000[,c(1,24)]
colnames(relp1000)[2] <- "P.value.1000"

diff <- merge(relp, relp1000, by = "Input.Sample")
diff <- diff[diff$P.value >=0.1 | diff$P.value.1000 >=0.1,]
