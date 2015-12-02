## Jinliang Yang
## July 17th, 2014
## plot the phenotypic traits


readPheno <- function(infile="akw_v2_GenSel_fullset.txt", outfile="akw_v2_GenSel_diallel.txt"){
  
  #setwd(pwd)
  pheno <- read.table(infile, header=TRUE)
  
  dpheno <- subset(pheno, pheno[,4] == "Diallel")
  message(sprintf("[ %s ] genotypes", nrow(dpheno)))
  ###output
  write.table(dpheno[, 1:3], outfile, sep="\t", row.names=FALSE, quote=FALSE)
  
  return(dpheno[, 1:3])
}


######
krn <- readPheno(infile="data/krn_v2_GenSel_fullset.txt", outfile="data/krn_v2_GenSel_diallel.txt")
cd <- readPheno(infile="data/cd_v2_GenSel_fullset.txt", outfile="data/cd_v2_GenSel_diallel.txt")
cl <- readPheno(infile="data/cl_v2_GenSel_fullset.txt", outfile="data/cl_v2_GenSel_diallel.txt")
cw <- readPheno(infile="data/cw_v2_GenSel_fullset.txt", outfile="data/cw_v2_GenSel_diallel.txt")
akw <- readPheno(infile="data/akw_v2_GenSel_fullset.txt", outfile="data/akw_v2_GenSel_diallel.txt")
tkw <- readPheno(infile="data/tkw_v2_GenSel_fullset.txt", outfile="data/tkw_v2_GenSel_diallel.txt")
kc <- readPheno(infile="data/kc_v2_GenSel_fullset.txt", outfile="data/kc_v2_GenSel_diallel.txt")


