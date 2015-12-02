# Jinliang Yang
# Purpose: control GenSel running
# date: 8.7.2011
# location: server.7


sh4Stepwise <- function(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
                      sh="CL_test.sh", email=TRUE,
                      geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                      pheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt"
                      ){
  #####
  setwd(pwd)
  inp <- gsub("sh", "inp", sh);
  log <- gsub("sh", "log", sh)
  
  Stepwise_inp(inp=inp, geno=geno, pheno=pheno, MaxR=0.8, MaxSNP=300, alpha=0.05)
  
  if(email){
    cat(paste("#Note:", Sys.time(), sep=" "),
        paste("GenSel4R", inp, ">", log),
        paste("python ~/bin/send_email.py -s ", "'", sh, "'", sep=""),
        file=sh, sep="\n", append=FALSE);
  }else{
    cat(paste("#Note:", Sys.time(), sep=" "),
        paste("GenSel4R", inp, ">", log),
        file=sh, sep="\n", append=FALSE);
  }
  
  
  
  message(paste("In this path: cd ", pwd, sep=""), "\n",
          paste("RUN: sh ", sh))
  
}

Stepwise_inp <- function(inp="CL_test.inp", 
                       geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                       pheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                       MaxR=0.8, MaxSNP=300, alpha=0.05){
  
    cat(paste("// gensel input file written", Sys.time(), sep=" "), 
        
        "analysisType StepWise",
        paste("inputMaxRsquared",  MaxR, sep=" "),
        paste("inputMaxMarkers",  MaxSNP, sep=" "),
        paste("alphaValue",  alpha, sep=" "),
        "",
        "// markerFileName",
        paste("markerFileName", geno, sep=" "), 
        "",
        "// phenotypeFileName",
        paste("phenotypeFileName", pheno, sep=" "),
        
        file=inp, sep="\n"
    )	
}
