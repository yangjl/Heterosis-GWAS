# Jinliang Yang
# Purpose: Prediction of SAM large and SAM small
# date: July.14.2014
# location: server.9


#####
source("lib/sh4GSCV.R")
setupCV_inp <- function(){
    ################ testrun for dom genotype
    traits <- c("krn", "cd", "cl", "cw", "akw", "tkw", "kc")
    genovar <- c(1.9, 7.6, 373, 91, 0.0017, 2515, 15485)
    residual <- c(0.3, 0.7, 105, 7, 0.0007, 178, 2258)
    
    for(modei in c("add", "dom")){
        for(ti in 1:7){ # traits
            for(j in 1:20){ ### change j to add more randomization 
                ####### dom and add
                for(k in 1:5){
                    mygeno <- paste("largedata/SNPdiallel/snp", traits[ti], 
                                    "_1m_", modei, "_gensel.newbin", sep="")
                    
                    inpfile <- paste0("largedata/scripts/", traits[ti],"_",modei,"_ran", j,"_cv", k, ".inp")
                    
                    mytrain <- paste("largedata/SNPdiallel/dcv1/pheno/", traits[ti],
                                     "_", j, "_train", k, ".txt", sep="")
                    mytest <- paste("largedata/SNPdiallel/dcv1/pheno/", traits[ti],
                                    "_", j, "_val", k, ".txt", sep="")
                    
                    GenSel_CV_inp(inp=inpfile, pi=0.9995, findsale ="no",
                                  geno=mygeno, 
                                  trainpheno=mytrain,
                                  testpheno=mytest,
                                  inmarker=NULL,
                                  chainLength=11000, burnin=1000, varGenotypic=genovar[ti], varResidual=residual[ti])
                }
            }
        }
    }
}
#######################
setupCV_inp()

files <- list.files(path="largedata/scripts", pattern="inp", full.names=TRUE)
logs <- gsub("inp", "log", files)
for(i in 1:length(files)){
    
    ### gerpIBD command goes into the sh
    sh1 <- paste("GenSel4R", files[i], ">", logs[i])
    cat(sh1, file=paste0("largedata/scripts/run_cv_job",i,".sh"), sep="\n", append=FALSE)
    
}

source("~/Documents/Github/zmSNPtools/Rcodes/set_arrayjob.R")
set_arrayjob(shid="largedata/run_array.sh",
             shcode=paste0("sh ", "largedata/scripts/run_cv_job$SLURM_ARRAY_TASK_ID.sh"),
             arrayjobs= "1-1400",
             wd=NULL, jobid="diallel", email="yangjl0930@gmail.com")
