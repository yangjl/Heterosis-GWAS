## Jinliang Yang
## July 17th, 2014
## plot the phenotypic traits


readPheno <- function(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
  infile="akw_v2_GenSel_fullset.txt", outfile="akw_v2_GenSel_diallel.txt"){
  
  setwd(pwd)
  pheno <- read.table(infile, header=TRUE)
  
  dpheno <- subset(pheno, pheno[,4] == "Diallel")
  message(sprintf("[ %s ] genotypes", nrow(dpheno)))
  ###output
  write.table(dpheno[, 1:3], outfile,
              sep="\t", row.names=FALSE, quote=FALSE)
  
  return(dpheno[, 1:3])
}


######
krn <- readPheno(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
                 infile="krn_v2_GenSel_fullset.txt", outfile="krn_v2_GenSel_diallel.txt")
cd <- readPheno(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
                 infile="cd_v2_GenSel_fullset.txt", outfile="cd_v2_GenSel_diallel.txt")
cl <- readPheno(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
                 infile="cl_v2_GenSel_fullset.txt", outfile="cl_v2_GenSel_diallel.txt")
cw <- readPheno(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
                 infile="cw_v2_GenSel_fullset.txt", outfile="cw_v2_GenSel_diallel.txt")
akw <- readPheno(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
                 infile="akw_v2_GenSel_fullset.txt", outfile="akw_v2_GenSel_diallel.txt")
tkw <- readPheno(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
                 infile="tkw_v2_GenSel_fullset.txt", outfile="tkw_v2_GenSel_diallel.txt")
kc <- readPheno(pwd="~/Documents/Heterosis_GWAS/HGWAS_proj/data/",
                 infile="kc_v2_GenSel_fullset.txt", outfile="kc_v2_GenSel_diallel.txt")


