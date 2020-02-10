pheno_score <- function (x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
{
  # setwd("/Users/qihua/Public/Phenotype-prediction")
  # source("pheno-score.R")
  all_pheno <-read.csv("selected_pheno.csv", stringsAsFactors = FALSE)
  ### need not adjust form
  csv <-read.csv("44_diy.csv", sep= "", stringsAsFactors = FALSE)
  phenofile <-readRDS("20191217_R2_age_sex.rds")


  predictors <- all_pheno[-1,]
  select_pheno <- all_pheno[1,]
  row.names(csv) <- csv$phenotype
  yinformat <- csv[ select_pheno, ]
  xinformat <-  csv[predictors, ]
  ##### obtain estimated beta and insample r2
  cor.x.y <- phenofile[select_pheno, predictors]
  cor.x <- phenofile[predictors, predictors]
  cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
  cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
  
  ebeta <- solve(cov.x) %*% cov.x.y
  # cov.y.hy <- cov.x.y %*% ebeta
  # var.hy <- t(sqrt(xinformat$var) * ebeta) %*% cor.x %*% (sqrt(xinformat$var) * ebeta)
  # cor.y.hy <- cov.y.hy^2/var.hy/yinformat$var

  ### used to identify the location of continuous variables
  con <- grep("_irnt", predictors)
  irnt_set <- predictors[con]
    
  value <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
  pheno_score <- sum(ebeta * value)
  return(pheno_score)
  
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

