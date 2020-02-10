selectvariable <- function (select_pheno)
{
  # setwd("/Users/qihua/Public/Prediction")
  # source("updated_variable_selection.R")
  library(stringr)
  library(Matrix)
  source("sojo.phenotype.function.R")
  csv <-read.csv("44_diy.csv", sep= "", stringsAsFactors = FALSE)
  phenofile <-readRDS("20191217_R2_age_sex.rds")
  
  phenotype.overlap <- intersect(csv$phenotype, colnames(phenofile))
  
  ##### obtain correlation matrix
  LD_phenomat <- phenofile[csv$phenotype, csv$phenotype]
  
  ##### select summary statistics concerning with and without ""

  yinformat <- csv[which(csv$phenotype == select_pheno), ]
  xinformat <-  csv[which(csv$phenotype != select_pheno), ]
  
  ##### obtain yxcorrelation and xxcorrelation
  vec <- which(colnames(LD_phenomat)== select_pheno)
  yxcor <- LD_phenomat[vec, -vec]
  xxcor <- LD_phenomat[-vec, -vec]
  df.sojo <- data.frame(var.x = xinformat$var, cor.x.y = yxcor)
  
  ##### When sum.stat.validation is not provided, LASSO result is provided #####
  res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = 10)
  # res$selected.markers <- c(res$selected.markers, "other")
  # reslist <- cbind(res$selected.markers, res$lambda.v)
  
  ##### When sum.stat.validation is same as sum.stat.discovery, insample R2 of LASSO result is provided #####
  res.insample <- sojo.phenotype(sum.stat.discovery = df.sojo, sum.stat.validation = df.sojo,
                                 cor.X = xxcor, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 10)
  insample_r2 <- res.insample$R2[2:11]
  
  reslist <- c( select_pheno, "sex", "age",res$selected.markers)
  # reslist <- c( select_pheno,res$selected.markers)
   write.csv(reslist, file = "selected_pheno.csv",row.names = FALSE)
   return(reslist)
}


