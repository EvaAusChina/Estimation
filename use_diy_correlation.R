# # setwd("/Users/qihua/Public/BS")
# library(Matrix)
# library(stringr)
# library("data.table")
# all_pheno <-read.csv("selected_pheno.csv", stringsAsFactors = FALSE)
# csv <-read.csv("44_diy.csv", sep= "", stringsAsFactors = FALSE)
# # phenofile <-readRDS("UNadjusted_R2.rds")
# 
# # predictors <- all_pheno[c(2,3),]
# # predictors <- all_pheno[-1,]
# # select_pheno <- all_pheno[1,]
# 
# 
# ##### input non White British data
# ts <- fread("/Users/qihua/Public/Prediction/5_new.tsv", header=TRUE, sep='\t', data.table=FALSE)
# 
# predictors <- all_pheno[-1,]
# select_pheno <- all_pheno[1,]
# select_data <- as.matrix(ts[, all_pheno[,1]])
# R2 <- cor(select_data,use = 'pairwise.complete.obs')
# yinformat <- csv[which(csv$phenotype == select_pheno), ]
# row.names(csv) <- csv$phenotype
# xinformat <-  csv[predictors, ]
# ##### obtain estimated beta and insample r2
# cor.x.y <- R2[select_pheno, predictors]
# cor.x <- R2[predictors, predictors]
# cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
# cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
# 
# ebeta <- solve(cov.x) %*% cov.x.y
# na_fix <- na.omit(select_data)
# p <- ncol(na_fix)
# q <- nrow(na_fix)
# score <- sweep(na_fix[,-1], 2, ebeta, "*")
# summary(lm(na_fix[,1]~ score[,1] ))
# lm_R2 <- c()
# for (i in 1:12) {
#   lm_summary <- summary(lm(na_fix[,1]~ score[,1:i] ))
#   lm_R2[i] <- lm_summary$r.squared
# }
# t <- data.frame(predictors, ebeta,lm_R2)
# 
# # > t
# # [[1]]
# # [1] 12
# # 
# # [[2]]
# # [1] 75066
# # 
# # [[3]]
# # [1] 0.04200363
# 
# predictors <- predictors[-c(1,2)]
# select_data <- as.matrix(ts[, all_pheno[-c(2,3),1]])
# R2 <- cor(select_data,use = 'pairwise.complete.obs')
# yinformat <- csv[which(csv$phenotype == select_pheno), ]
# row.names(csv) <- csv$phenotype
# xinformat <-  csv[predictors, ]
# ##### obtain estimated beta and insample r2
# cor.x.y <- R2[select_pheno, predictors]
# cor.x <- R2[predictors, predictors]
# cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
# cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
# ebeta <- solve(cov.x) %*% cov.x.y
# 
# lm_R2 <- c()
# na_fix <- na.omit(select_data)
# p <- ncol(na_fix)
# q <- nrow(na_fix)
# score <- sweep(na_fix[,-1], 2, ebeta, "*")
# for (i in 1:length(predictors)) {
#   lm_summary <- summary(lm(na_fix[,1]~ score[,1:i] ))
#   lm_R2[i] <- lm_summary$r.squared
# }
# t <- data.frame(predictors, ebeta,lm_R2)
# 
# 
# 
### correlation
indi <- fread("44_new.tsv", header=TRUE, sep='\t', data.table=FALSE)
ts <- fread("5_new.tsv", header=TRUE, sep='\t', data.table=FALSE)
all_pheno <-read.csv("selected_pheno.csv", stringsAsFactors = FALSE)

# indi <- indi[1:50000,]
library(stringr)
library(Matrix)
source("sojo.phenotype.function.R")

csvv <-read.csv("44_diy.csv", sep= "", stringsAsFactors = FALSE)
row.names(csvv) <- csvv$phenotype
csv <- csvv[which(csvv$phenotype != c("3526_irnt")), ]
# csv2 <- csv1[which(csv1$phenotype != c("23111_irnt")), ]
# csv <- csv2[which(csv2$phenotype != c("23099_irnt")), ]


a <-readRDS("20191217_R2_age_sex.rds")
phenofile <- a[csv$phenotype,csv$phenotype]
  ##### obtain yxcorrelation and xxcorrelation
  vec <- which(colnames(phenofile)== select_pheno)
  yxcor <- phenofile[-vec, vec]
  xxcor <- phenofile[-vec, -vec]
  xinformat <- csv[-vec,]
  yinformat <- csv[vec, ]
  df.sojo <- data.frame(var.x = xinformat$var, cor.x.y = yxcor)
  
  res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = 10)
  res.insample <- sojo.phenotype(sum.stat.discovery = df.sojo, sum.stat.validation = df.sojo,
                                 cor.X = xxcor, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 11)
  
   select_pheno <- "1807_irnt"
  select_data <- as.matrix(ts[, c(select_pheno, "sex", "age",res$selected.markers)])

    ### with 17 variables,15 selected, 
  R2 <- cor(select_data,use = 'pairwise.complete.obs')
  
    
#### predictors are selected phenos(<= 12) 
  {  
pheno_log <- function(num){
  res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = num)
  # predictors <- c("sex", "age",res$selected.markers)
  predictors <- res$selected.markers
  cor.x.y <- R2[select_pheno, predictors]
  cor.x <- R2[predictors, predictors]
  xinformat <-  informat[predictors, ]
  yinformat <- informat[select_pheno, ]
  if(num == 1){
    cov.x <- sqrt(xinformat$var) * cor.x * sqrt(xinformat$var)
    cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
    be <- solve(cov.x) %*% cov.x.y
    cov.y.hy2 <- cov.x.y * be
    var.hy2 <- (sqrt(xinformat$var) * be)* cor.x*(sqrt(xinformat$var) * be)
    cor.y.hy2 <- cov.y.hy2^2/var.hy2/yinformat$var
  }else{
  cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
  cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
  be <- solve(cov.x) %*% cov.x.y
  cov.y.hy2 <- cov.x.y %*% be
  var.hy2 <- t(sqrt(xinformat$var) * be) %*% cor.x %*% (sqrt(xinformat$var) * be)
  cor.y.hy2 <- cov.y.hy2^2/var.hy2/yinformat$var}
  ## out_r2
  lm_R2 <- c()
 
  # va <-  c("1807_irnt", "sex", "age",res$selected.markers)
  va <-  c(select_pheno, res$selected.markers)
  b <- be
  ta <- as.matrix(indi[, va])
  na_fix <- na.omit(ta)
  p <- ncol(na_fix)
  q <- nrow(na_fix)
 
  if (num == 1){
    score <- na_fix[,-1] * as.vector(b)
    lm_summary <- summary(lm(na_fix[,1]~ score ))
    lm_R2[1] <- lm_summary$r.squared
  }
  else{
    score <-  sweep(na_fix[,-1], 2, b, "*")
    for (i in 2:length(va[-1])) {
      lm_summary <- summary(lm(na_fix[,1]~ score[,1:i] ))
      lm_R2[i] <- lm_summary$r.squared
    }
  }
  t <- list(va[length(va)], p-1, q,lm_R2[length(lm_R2)], be[length(be)],cor.y.hy2)
  return(t)
}

pheno <- c()
var_num <- c()
sam_size <- c()
r2 <- c()
coeff <- c()
in_sam <-c()
for (i in 1:10) {
  m <- pheno_log(i)
  pheno[i] <- m[[1]]
  var_num[i] <- m[[2]]
  sam_size[i] <- m[[3]]
  r2[i] <- m[[4]]
  coeff[i] <- m[[5]]
  in_sam[i] <- m[[6]]
}

s <- data.frame(pheno, var_num, sam_size, r2,coeff,in_sam,stringsAsFactors = FALSE)
colnames(s) <- c("pheno","var_num", "sam_size", "out-of-sample-r2","coefficient","in-sample-r2")

# write.table(s, file = "out_sample_without 6177_100.csv",row.names = FALSE)
  }
  
  
### predictors are selected phenos + sex
{
predictors <- c("sex", res$selected.markers)
cor.x.y <- R2[select_pheno, predictors]
cor.x <- R2[predictors, predictors]
xinformat <-  csv
row.names(xinformat) <- xinformat$phenotype
xinformat <-  xinformat[predictors, ]
  cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
  cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
  be <- solve(cov.x) %*% cov.x.y
  va <-  c(select_pheno, "sex",res$selected.markers)
  b <- be
  mark <- NULL
  select_data <- as.matrix(indi[, va])
  na_fix <- na.omit(select_data)
  p <- ncol(na_fix)
  q <- nrow(na_fix)
    score <- sweep(na_fix[,-1], 2, b, "*")
      lm_summary <- summary(lm(na_fix[,1]~ score[,1:11] ))
      lm_R2 <- lm_summary$r.squared
      sex_info <- c("sex", 11, q, lm_R2,b)
}      

### predictors are selected phenos + sex + age
{
      predictors <- c("sex", "age",res$selected.markers)
      cor.x.y <- R2[select_pheno, predictors]
      cor.x <- R2[predictors, predictors]
      xinformat <-  csv
      row.names(xinformat) <- xinformat$phenotype
      xinformat <-  xinformat[predictors, ]
      cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
      cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
      be <- solve(cov.x) %*% cov.x.y
      va <-  c(select_pheno, "sex", "age",res$selected.markers)
      b <- be
      mark <- NULL
      select_data <- as.matrix(indi[, va])
      na_fix <- na.omit(select_data)
      p <- ncol(na_fix)
      q <- nrow(na_fix)
      score <- sweep(na_fix[,-1], 2, b, "*")
      lm_summary <- summary(lm(na_fix[,1]~ score[,1:12] ))
      lm_R2 <- lm_summary$r.squared
      sa_info <- c("age", 12, q, lm_R2,b[2])
}

### predictor is only sex
{
lm_R2 <- c()
cor.x.y <- R2[select_pheno, "sex"]
cor.x <- R2["sex", "sex"]
xinformat <-  csv
row.names(xinformat) <- xinformat$phenotype
xinformat <-  xinformat["sex", ]
cov.x <- sqrt(xinformat$var) * cor.x * sqrt(xinformat$var)
cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
be <- solve(cov.x) %*% cov.x.y
sex_set <- as.matrix(ts[,c(select_pheno, "sex")])
na_fix <- na.omit(sex_set)
p <- ncol(na_fix)
q <- nrow(na_fix)
score <- na_fix[,-1] * as.vector(be)
lm_summary <- summary(lm(na_fix[,1]~ score ))
lm_R2 <- lm_summary$r.squared
sex_info <- c("sex", 1, q, lm_R2)
} 

### predictors are sex + age
{
lm_R2 <- c()
sa <- c("sex","age")
sa_set <- as.matrix(ts[,c(select_pheno,"sex","age")])
cor.x.y <- R2[select_pheno, sa]
cor.x <- R2[sa, sa]
xinformat <-  csv
row.names(xinformat) <- xinformat$phenotype
xinformat <-  xinformat[sa, ]
cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
be <- solve(cov.x) %*% cov.x.y

na_fix <- na.omit(sa_set)
p <- ncol(na_fix)
q <- nrow(na_fix)
score <- sweep(na_fix[,-1], 2, be, "*")
  lm_summary <- summary(lm(na_fix[,1]~ score[,1:2] ))
  lm_R2 <- lm_summary$r.squared
sa_info <- c("age", 2, q, lm_R2)
}

  all_in <- rbind.data.frame(s,sex_info, sa_info)
 write.table(all_in, file = "outsample_R2_3.csv",row.names = FALSE)
  