# Jinliang Yang
# Purpose: control GenSel running
# date: Feb.10.2016
# location: server.9

############################
GenSel_CV_inp <- function(inp="CL_test.inp", pi=0.995, findsale ="no",
                       geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                       trainpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                       testpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                       inmarker="/Users/yangjl/Documents/linkage",
                       chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2){
    
    cat(paste("// gensel input file written", Sys.time(), sep=" "), 
        
        "analysisType Bayes",
        "bayesType BayesC",
        paste("chainLength", chainLength, sep=" "),
        paste("burnin", burnin=burnin, sep=" "),
        paste("probFixed", pi, sep=" "),
        
        paste("varGenotypic",  varGenotypic, sep=" "),
        paste("varResidual",  varResidual, sep=" "),
        "nuRes 10",
        "degreesFreedomEffectVar 4",
        "outputFreq 100",
        "seed 1234",
        "mcmcSamples yes",
        "plotPosteriors no",
        paste("FindScale", findsale),
        "modelSequence no",
        "isCategorical no",
        
        "",
        "// markerFileName",
        paste("markerFileName", geno, sep=" "), 
        "",
        "// train phenotypeFileName",
        paste("trainPhenotypeFileName", trainpheno, sep=" "),
        "",
        "// test phenotypeFileName",
        paste("testPhenotypeFileName", testpheno, sep=" "),
        
        "// includeFileName",
        #if(!is.null(inmarker)){
         #   paste("includeFileName", inmarker, sep=" ")
        #}
       
        file=inp, sep="\n"
    )	
}
