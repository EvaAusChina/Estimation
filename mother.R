
indi <- fread("5_new.tsv", header=TRUE, sep='\t', data.table=FALSE)
ts <- fread("44_new.tsv", header=TRUE, sep='\t', data.table=FALSE)
all_pheno <-read.csv("selected_pheno.csv", stringsAsFactors = FALSE)

mo <- indi[,"3526_irnt"]
max(mo)

library(stringr)
library(Matrix)
source("sojo.phenotype.function.R")
csvv <-read.csv("44_diy.csv", sep= "", stringsAsFactors = FALSE)
row.names(csvv) <- csvv$phenotype
csv <- csvv[which(csvv$phenotype != c("1807_irnt")), ]
# csv2 <- csv1[which(csv1$phenotype != c("23111_irnt")), ]
# csv <- csv2[which(csv2$phenotype != c("23099_irnt")), ]
select_pheno <- "3526_irnt"


a <-readRDS("20191217_R2_age_sex.rds")
phenofile <- a[csv$phenotype,csv$phenotype]
##### obtain yxcorrelation and xxcorrelation
vec <- which(colnames(phenofile)== select_pheno)
yxcor <- phenofile[-vec, vec]
xxcor <- phenofile[-vec, -vec]
xinformat <- csv[-vec,]
yinformat <- csv[vec, ]
df.sojo <- data.frame(var.x = xinformat$var, cor.x.y = yxcor)

res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = 14)
select_data <- as.matrix(ts[, c(select_pheno, "sex", "age",res$selected.markers)])

### with 17 variables,15 selected, 
R2 <- cor(select_data,use = 'pairwise.complete.obs')


#### predictors are selected phenos(<= 12) 
{  
  pheno_log <- function(num){
    res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = num)
    # predictors <- c("sex", "age",res$selected.markers)
    predictors <- c(res$selected.markers)
    cor.x.y <- phenofile[select_pheno, predictors]
    cor.x <- phenofile[predictors, predictors]
    xinformat <-  csv
    row.names(xinformat) <- xinformat$phenotype
    xinformat <-  xinformat[predictors, ]
    if(num == 1){
      cov.x <- sqrt(xinformat$var) * cor.x * sqrt(xinformat$var)
      cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
      be <- solve(cov.x) %*% cov.x.y
    }else{
      cov.x <- diag(sqrt(xinformat$var)) %*% cor.x %*% diag(sqrt(xinformat$var))
      cov.x.y <- cor.x.y * sqrt(yinformat$var) * sqrt(xinformat$var)
      be <- solve(cov.x) %*% cov.x.y}
    ## out_r2
    lm_R2 <- c()
    # for (i in 2:length(predictors)) {
    #   if(predictors[i] %in% colnames(ts) == FALSE){
    #     predictors[i] <- NA
    #     be[i-1] <- NA
    #   }
    # }
    # 
    # va <- na.omit(predictors)
    # b <- na.omit(be)
    # mark <- as.vector(na.action(b))
    # va <-  c("1807_irnt", "sex", "age",res$selected.markers)
    va <-  c(select_pheno, res$selected.markers)
    b <- be
    mark <- NULL
    select_data <- as.matrix(indi[, va])
    na_fix <- na.omit(select_data)
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
    t <- list(va[length(va)], p-1, q,lm_R2[length(lm_R2)], mark, be[length(be)])
    return(t)
  }
  
  pheno <- c()
  var_num <- c()
  sam_size <- c()
  r2 <- c()
  coeff <- c()
  for (i in 1:14) {
    m <- pheno_log(i)
    pheno[i] <- m[[1]]
    var_num[i] <- m[[2]]
    sam_size[i] <- m[[3]]
    r2[i] <- m[[4]]
    om <- m[[5]]
    coeff[i] <- m[[6]]
  }
  if(class(om) != "numeric"){
    s <- data.frame(pheno, var_num, sam_size, r2,coeff,stringsAsFactors = FALSE)
  }else{
    s <- data.frame(pheno[-om], var_num[-om], sam_size[-om], r2[-om],coeff[-om])
  }
  colnames(s) <- c("pheno","var_num", "sam_size", "r2","coefficient")
  # write.table(s, file = "out_sample_without 6177_100.csv",row.names = FALSE)
}