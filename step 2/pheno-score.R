pheno_score <- function (x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
{
# setwd("/Users/qihua/Public/Phenotype-prediction")
# source("pheno-score.R")
all_pheno <-read.csv("selected_pheno.csv", stringsAsFactors = FALSE)
pheno_summary <-read.csv("summary_pheno.csv", stringsAsFactors = FALSE)
csv <-read.csv("updated_phenofile.csv", stringsAsFactors = FALSE)
  
phenofile <-readRDS("20190901_R2.rds")
predictors <- all_pheno[-1,]
select_pheno <- all_pheno[1,]

yinformat <- pheno_summary[which(pheno_summary$phenotype == select_pheno), ]
xinformat <-  pheno_summary[match(predictors, pheno_summary$phenotype), ]
LD_phenomat <- phenofile[predictors, predictors]

##### obtain estimated beta and insample r2
cor.x.y <- phenofile[select_pheno, predictors]
cor.x <- LD_phenomat
cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)

ebeta <- solve(cov.x) %*% cov.x.y
# cov.y.hy <- cov.x.y %*% ebeta
# var.hy <- t(sqrt(xinformat$var) * ebeta) %*% cor.x %*% (sqrt(xinformat$var) * ebeta)
# cor.y.hy <- cov.y.hy^2/var.hy/yinformat$var

### handle irnt data
# load("bd.ukb11243.RData")
# loca <- which(colnames(bd) == "f.21000.0.0")
# trybd <- bd[which(bd[,loca] == "1001"),]
# loc <- which(colnames(trybd) == "f.22006.0.0")
# bd1 <- trybd[,-loc]
# rename <- bd1[, -1]
# re1 <- gsub("^f.", "", colnames(rename))
# re2 <- gsub(".[0-9]+.[0-9]+$", "", re1)
# colnames(rename) <- re2
# colnames(rename) <- paste(colnames(rename), "_irnt", sep = "")
# save(rename,file="bd_irnt.Rdata")
load("bd_irnt.Rdata")

irnt <- function(cts_variable) {
  set.seed(1234) 
  n_cts <- length(which(!is.na(cts_variable)))
  quantile_cts <- (rank(cts_variable, na.last = "keep", ties.method = "random") - 0.5) / n_cts
  cts_IRNT <- qnorm(quantile_cts)	
  return(cts_IRNT)
}
vice_irnt <- function(cts_IRNT) {
  set.seed(1234)
  quantile_cts <- pnorm(cts_IRNT)
  n_cts <- length(which(!is.na(cts_IRNT)))
  rank <- quantile_cts * n_cts 
  rit <- c(rank[1], n_cts, quantile_cts * 100)
  return(rit)
}

value <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
for (i in 1:length(predictors)) {
  if(csv[which(csv$phenotype == predictors[i]), 3] == "continuous_irnt"){
    if(predictors[i] %in% colnames(rename)){
    irnt_set <- c(value[i], rename[,predictors[i]])
     ir <- irnt(irnt_set)
     value[i] <- ir[i]
    }else{
      value[i] <- NA
      ebeta[i] <- NA
    }
    }
}

pheno_score <- sum(na.omit(ebeta * value))

if(csv[which(csv$phenotype == select_pheno), 3] == "binary"){
  cat("you have", pheno_score, "probability of the phenotype")
}else{
  con <- read.csv("continuous_set.csv", sep = "", stringsAsFactors = FALSE)
  colnames(con) <- gsub("X", "", colnames(con))
  virnt_set <- c(pheno_score, con[, select_pheno])
  vi <- vice_irnt(virnt_set)
  cat("You rank", vi[1], "in", vi[2], "people.","\n", "That is", vi[3], "%","\n")
 }
}
