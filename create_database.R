### 5万估计(white)，5万验证
indi <- fread("44_new.tsv", header=TRUE, sep='\t', data.table=FALSE)
ts <- fread("5_new.tsv", header=TRUE, sep='\t', data.table=FALSE)
# ran <- 1:57160
# qa <- sample(x = ran, size = 57160, replace = TRUE)
# indi <- indi[qa,]
# ts <- ts[qa,]
# lap_set <- ts
# meanfun <- function(x){
#   a <- na.omit(x)
#   return(mean(a))
# }
# varfun <- function(x){
#   a <- na.omit(x)
#   return(sd(a))
# }
# me <- apply(lap_set, 2, meanfun)
# va <- apply(lap_set, 2, varfun)
# mv <- data.frame(colnames(lap_set), me, va)
# colnames(mv) <- c("phenotype", "mean", "var")
# csvv <- mv

csvv <-read.csv("44_diy.csv", sep= "", stringsAsFactors = FALSE)
row.names(csvv) <- csvv$phenotype
csv <- csvv[which(csvv$phenotype != c("1807_irnt")), ]
informat <-  csv
row.names(informat) <- informat$phenotype
xinformat <-  informat[predictors, ]
yinformat <- informat[select_pheno, ]

select_data <- as.matrix(ts[, c(select_pheno,res$selected.markers)])
R2 <- cor(select_data,use = 'pairwise.complete.obs')
predictors <- res$selected.markers
cor.x.y <- R2[select_pheno, predictors]
cor.x <- R2[predictors, predictors]

yxcor <- cor.x.y
xxcor <- cor.x
df.sojo <- data.frame(var.x = xinformat$var, cor.x.y = yxcor)
res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = 10)
res.insample <- sojo.phenotype(sum.stat.discovery = df.sojo, sum.stat.validation = df.sojo,
                               cor.X = xxcor, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 10)




### 440000 white
train <- ts[,c(select_pheno,res$selected.markers)] 
valid <- indi[,c(select_pheno,res$selected.markers)]

simu <- function(){
si_cov2 <- cov(train,use = 'pairwise.complete.obs')
si_cor2 <- cor(train,use = 'pairwise.complete.obs')
me2 <- colMeans(train, na.rm = TRUE)
si_tr <- rmvnorm(nrow(train), me2, si_cov2)
dim(si_tr)
### 50000 non-white
si_cov <- cov(valid,use = 'pairwise.complete.obs')
si_cor <- cor(valid,use = 'pairwise.complete.obs')
me1 <- colMeans(valid, na.rm = TRUE)
si_va <- rmvnorm(nrow(indi), me1, si_cov)
dim(si_va)
colnames(si_tr) <- colnames(train)
colnames(si_va) <- colnames(train)

### input
R2 <- si_cor2
meanfun <- function(x){
  a <- na.omit(x)
  return(mean(a))
}
varfun <- function(x){
  a <- na.omit(x)
  return(sd(a))
}
length(me)
me <- apply(si_tr, 2, meanfun)
va <- apply(si_tr, 2, varfun)
mv <- data.frame(colnames(si_tr), me, va)
colnames(mv) <- c("phenotype", "mean", "var")
informat <- mv
row.names(informat) <- informat$phenotype
indi <- si_va

yxcor <- R2[select_pheno, res$selected.markers]
xxcor <- R2[res$selected.markers, res$selected.markers]
xinformat <-  informat[res$selected.markers, ]
yinformat <- informat[select_pheno, ]
df.sojo <- data.frame(var.x = xinformat$var, cor.x.y = yxcor)

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
    cor.y.hy2 <- cov.y.hy2^2/var.hy2/yinformat$var
    }
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
# s <- data.frame(pheno, var_num, sam_size, r2,coeff,in_sam,stringsAsFactors = FALSE)
# colnames(s) <- c("pheno","var_num", "sam_size", "out-of-sample-r2","coefficient","in-sample-r2")
return(r2)
}

  library(plotly)
p<-plot_ly(y = simu(), x= 1:10 , type="scatter", mode="markers+lines", name = "1")

for (i in 2:20) {
p<-add_trace(p, y = simu(), x=1:10, type="scatter", mode="markers+lines", name = i)
}









