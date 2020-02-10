# source("/Users/qihua/Public/Phenotype-prediction/creat_new_phenofile_csv.R")
library(Matrix)
library(stringr)
csv <- read.csv("/Users/qihua/Public/phenotypes.both_sexes.csv", sep = ";", stringsAsFactors = FALSE)
bin_info1 <- csv[csv[,3]=="binary",c(1,2,3,5,6,7,8)]
bin_info <- bin_info1[!(str_detect(bin_info1[,2], "Illnesses")),]
irnt_info <- csv[csv[,3]=="continuous_irnt",c(1,2,3,5,6,7,8)]
csv <- rbind(bin_info, irnt_info)
cs <- csv[csv[,3]=="binary",]
csvdata <- cs[, c(1,4,5,6,7)]
ncase <- as.numeric(csvdata[,5])
ncontrol <- as.numeric(csvdata[,4])
n1 <- as.numeric(csvdata[,2])
##### obtain mean and var of binary data
nmean <- ncase / n1
nvar <- ncontrol * ncase / (n1 * n1)
csvdata[,6] <- nmean
csvdata[,7] <- nvar
colnames(csvdata)[6] <- "mean"
colnames(csvdata)[7] <- "var"
csvb <- csvdata[,c(1,6,7)]

##### obtain mean and var of continuous data
csvc <- csv[csv[,3]=="continuous_irnt",]
len <- nrow(csvc)
csvc[,2] <- rep(0, len)
csvc[,3] <- rep(1, len)
csvc <- csvc[,c(1,2,3)]
colnames(csvc) <- colnames(csvb)
csvn <- rbind(csvb, csvc)

tsv_data <- fread("WAS/results/output..tsv", header=TRUE, sep='\t', data.table=FALSE)
pheno_data <- tsv_data[,-c(1,2,3)]
judge_set <- c()
binary_set <- pheno_data[,which(judge_set)]
binary_set1 <- ifelse(binary_set == TRUE, 1, 0)
for (i in 1:length(colnames(pheno_data))) {
  if(class(pheno_data[,i]) == "logical"){
    judge_set[i] = TRUE
  }else{
    judge_set[i] = FALSE
  }}
continuous_set <- pheno_data[,-which(judge_set)]
colnames(continuous_set) <- paste(colnames(continuous_set), "_irnt", sep = "")

total_info <- csvn[which(csvn$phenotype %in% colnames(pheno_data)),  ]

write.csv(total_info, file = "/Users/qihua/Public/Phenotype-prediction/updated_phenofile.csv",row.names = FALSE)

### unfinished

