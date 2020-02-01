selectvariable <- function (select_pheno)
{
  library(stringr)
  library(Matrix)
  source("sojo.phenotype.function.R")
  csv <-read.csv("updated_phenofile.csv", stringsAsFactors = FALSE)
  csvn <- csv[csv[,3]=="binary",]
  csvdata <- csvn[, c(1,4,5,6,7)]
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
  
  ##### seek overlap between correlation matrix and phenotype table
  phenofile <-readRDS("G:/phenotype estimation/sojo/20190901_R2.rds")
  phenotype.overlap <- intersect(csvn$phenotype, colnames(phenofile))
  colnames(phenofile)
  
  ##### obtain correlation matrix
  LD_phenomat <- phenofile[phenotype.overlap, phenotype.overlap]
  
  ##### select summary statistics concerning with and without ""
  informat <-  csvn[csvn$phenotype %in% phenotype.overlap,]
  yinformat <- subset(informat, phenotype == select_pheno, select = c(phenotype, mean, var))
  xinformat <-  subset(informat, phenotype != select_pheno, select = c(phenotype, mean, var))
  
  ##### obtain yxcorrelation and xxcorrelation
  vec <- which(colnames(LD_phenomat)== select_pheno)
  yxcor <- LD_phenomat[colnames(LD_phenomat)== select_pheno, -vec]
  xxcor <- LD_phenomat[-vec, -vec]
  df.sojo <- data.frame(var.x = xinformat$var, cor.x.y = yxcor)
  
  ##### When sum.stat.validation is not provided, LASSO result is provided #####
  res <- sojo.phenotype(sum.stat.discovery = df.sojo, cor.X = xxcor, v.y = yinformat$var, nvar = 10)
  # res$selected.markers <- c(res$selected.markers, "other")
  # reslist <- cbind(res$selected.markers, res$lambda.v)
  
  ##### When sum.stat.validation is same as sum.stat.discovery, insample R2 of LASSO result is provided #####
  res.insample <- sojo.phenotype(sum.stat.discovery = df.sojo, sum.stat.validation = df.sojo,
                                 cor.X = xxcor, v.y = yinformat$var, v.y.validation = yinformat$var, nvar = 10)
  
  insample_r2 <- res.insample$R2[2:21]
  
  reslist <- cbind(res$selected.markers, insample_r2)
  colnames(reslist) <- c("variable", "r2")
   print(reslist)
}


