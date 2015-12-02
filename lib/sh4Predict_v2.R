### Jinliang Yang
### generate shell file for Bayes Predict

sh4Predict_v2 <- function(pwd="~/Documents/Heterosis_GWAS/Method/GenSel/testrun", 
                       sh="CL_test.sh", inp="some.inp",
                       geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                       pheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                       map="/Users/yangjl/Documents/linkage.map",
                       mrkRes="", shappend=FALSE){
  #####
  setwd(pwd)
  log <- gsub("inp", "log", inp)
  Predict_inp(inp=inp, geno=geno, pheno=pheno, map=map, mrkRes=mrkRes)
  
  cat(paste("#Note:", Sys.time(), sep=" "),
      paste("GenSel4R", inp, ">", log),
      paste("#python ~/bin/send_email.py -s ", "'", sh, "'", sep=""),
      file=sh, sep="\n", append=shappend);
  
  message(sprintf("In this path: [ %s ],\nRUN: [ sh %s ]", pwd, sh))
  
}


Predict_inp <- function(inp="",
                        geno="",
                        pheno="",
                        map="",
                        mrkRes=""){
  cat(paste("// gensel input file written", Sys.time(), sep=" "), 
      
      "analysisType Predict",
      "windowWidth 5",
      "#linkageMap AGPv2",
      "#addMapInfoToMarkers yes",
      
      "// markerFileName",
      paste("markerFileName", geno, sep=" "), 
      "",
      "// phenotypeFileName",
      paste("phenotypeFileName", pheno, sep=" "),
      "",
      "// mapOrderFileName",
      paste("// mapOrderFileName", map, sep=" "),
      "",
      "// markerSolFileName",
      paste("markerSolFileName", mrkRes, sep=" "),
      
      file=inp, sep="\n")
}