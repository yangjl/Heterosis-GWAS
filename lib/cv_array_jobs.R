### Jinliang Yang








source("~/Documents/Github/zmSNPtools/Rcodes/set_arrayjob.R")
setup_gerpibd_array_7traits <- function(
  outdir="slurm-scripts/cv_b0/", jobbase="run_cv_job", jobid =1,
  kfile_path="largedata/snpeff/BPH/",
  genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0"){
  #jobid: the start number of jobid
  
  traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
  dir.create(outdir, showWarnings = FALSE)
  jobstart <- jobid
  for(j in 1:7){
    shid <- paste0(outdir, "/", jobbase, jobid, ".sh")
    ### gerpIBD command goes into the sh
    sh1 <- paste0("gerpIBD -d largedata/IBD/allsnps_11m_IBD.bed -s largedata/SNP/allsnps_11m.dsf5 ",
                  "-g ", genobase, ".csv -f ", kfile_path, "/", traits[j], "_k.txt ",
                  "-o ", genobase, "_", traits[j],
                  " -t k")
    cat(sh1, file=shid, sep="\n", append=FALSE)
    message(sprintf("###>>> traits [ %s ]; total jobs [ %s ]", traits[j], jobid))
    jobid <- jobid+1
  }
  jobend <- jobid -1
  message(sprintf("###>>> setup array jobs: [ %s - %s]", jobstart, jobend))
  set_arrayjob(shid=paste0(outdir, "/", jobbase, ".sh"),
               shcode=paste0("sh ", outdir, "/", jobbase, "$SLURM_ARRAY_TASK_ID.sh"),
               arrayjobs= paste0("1-", jobend),
               wd=NULL, jobid=jobbase, email="yangjl0930@gmail.com")
  
}

###########################################
source("lib/slurm4gerpIBDCV.R")
setup_newbin_array <- function(
  ### note: it is for 7 traits with 3 modes for one random shuffling or real data
  genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0",
  phenobase="largedata/pheno/CV5fold",
  jobdir="slurm-scripts/get_newbin", inpbase= "cs0",
  jobbase="run_newbin_job", jobid =1){
  
  ### prior information
  wd <- getwd()
  #test run of the 66 diallel of trait per se with additive model
  ti <- tolower(c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW"))
  res <- c(0.38, 0.46, 0.46, 15, 88, 41, 0.64)
  gen <- c(0.18, 5.1, 6.0, 123, 65, 377, 0.82)
  
  dir.create(jobdir, showWarnings = FALSE)
  shcommand <- c()
  for(myti in 1:7){
    for(mode in c("a2", "d2", "h2")){
      ### the first one use gs
      cv <- 1
      sp <- 1
      myinp <- paste0(jobdir, "/", inpbase, "_", ti[myti], "_", mode,"_cv",cv, "_sp",sp, ".inp")
      GS_cv_inp(
        inp= myinp, pi=0.995,
        # largedata/SNP/geno_b0_cs/gerpIBD_b0_cs", i, "_", traits[j]
        geno= paste0(wd, "/", genobase, "_", ti[myti], "_", mode, ".gs"), 
        
        #phenobase
        trainpheno= paste0(wd, "/", phenobase, "/", ti[myti], "_train", cv, "_sp", sp, ".txt"),
        testpheno= paste0(wd, "/", phenobase, "/", ti[myti], "_test", cv, "_sp", sp, ".txt"),
        chainLength=100, burnin=10, varGenotypic=gen[myti], varResidual=res[myti]
      )
      shcommand <- c(shcommand, paste("GenSel4R", myinp))
    }
  }
  #################
  jobstart = jobid
  for(i in 1:length(shcommand)){
    cat(shcommand[i], file=paste0(jobdir, "/", jobbase, jobid, ".sh"), sep="\n", append=FALSE)
    jobid <- jobid + 1
  }
  jobend <- jobid -1
  message(sprintf("###>>> setup array jobs: [ %s - %s]", jobstart, jobend))
  set_arrayjob(shid=paste0(jobdir, "/", jobbase, ".sh"),
               shcode=paste0("sh ", jobdir, "/", jobbase, "$SLURM_ARRAY_TASK_ID.sh"),
               arrayjobs= paste0("1-", jobend),
               wd=NULL, jobid=jobbase, email="yangjl0930@gmail.com")
  
}
#newbin_array_7traits_3modes(genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0",
#                            jobdir="slurm-scripts/get_newbin/", jobbase="run_newbin_job")

  
  

setup_gensel_array <- function(outdir="slurm-scripts/cv_b2", jobbase="run_gs_job", jobid=1,
                               inpbase="slurm-scripts/cv_b2/cs0", 
                               phenobase="largedata/pheno/CV5fold",
                       genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0"){
  
  traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
  dir.create(outdir, showWarnings = FALSE)
  
  jobstart <- jobid
  for(sp in 1:100){# 10 sub-sampling
    
    sh1 <- c() #=> 105 shell commands
    for(i in 1:7){# 7 traits
      mysh <- cv15_mode_cv(inpbase=inpbase, genobase=genobase, phenobase=phenobase, myti=i, sp)
      sh1 <- c(sh1, mysh)
    }
    
    sh2 <- c(paste0("rm ", inpbase, "_", "*sp", sp, ".mcmcSamples1"),
             paste0("rm ", inpbase, "_", "*sp", sp, ".mrkRes1"),
             paste0("rm ", inpbase, "_", "*sp", sp, ".inp"),
             paste0("rm ", inpbase, "_", "*sp", sp, ".out1"),
             paste0("rm ", inpbase, "_", "*sp", sp, ".cgrRes1"))
    
    shid <- paste0(outdir, "/", jobbase, jobid, ".sh")
    cat(c(sh1, sh2), file=shid, sep="\n", append=FALSE)
    jobid <- jobid +1
  }
  jobend <- jobid - 1
  
  #################
  message(sprintf("###>>> setup gensel array jobs: [ %s - %s]", jobstart, jobend))
  set_arrayjob(shid=paste0(outdir, "/", jobbase, ".sh"),
               shcode=paste0("sh ", outdir, "/", jobbase, "$SLURM_ARRAY_TASK_ID.sh"),
               arrayjobs= paste0("1-", jobend),
               wd=NULL, jobid=jobbase, email="yangjl0930@gmail.com")
}



############
cv15_mode_cv <- function(inpbase="slurm-scripts/cv_b0/cs0", 
                         genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0", 
                         phenobase="largedata/pheno/CV5fold",
                         myti, sp){
  ## myti: ith trait; sp: nth sub sampling
  
  ### prior information
  wd <- getwd()
  #test run of the 66 diallel of trait per se with additive model
  ti <- tolower(c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW"))
  res <- c(0.38, 0.46, 0.46, 15, 88, 41, 0.64)
  gen <- c(0.18, 5.1, 6.0, 123, 65, 377, 0.82)
  
  shcommand <- c()
  for(mode in c("a2", "d2", "h2")){
    for(cv in 1:5){
      myinp <- paste0(inpbase, "_", ti[myti], "_", mode,"_cv",cv, "_sp",sp, ".inp")
      GS_cv_inp(
        inp= myinp, pi=0.995,
        # largedata/SNP/geno_b0_cs/gerpIBD_b0_cs", i, "_", traits[j]
        geno= paste0(wd, "/", genobase, "_", ti[myti], "_", mode, ".gs.newbin"), 
        trainpheno= paste0(wd, "/", phenobase, "/", ti[myti], "_train", cv, "_sp", sp, ".txt"),
        testpheno= paste0(wd, "/", phenobase, "/", ti[myti], "_test", cv, "_sp", sp, ".txt"),
        chainLength=11000, burnin=1000, varGenotypic=gen[myti], varResidual=res[myti])
        
        shcommand <- c(shcommand, paste("GenSel4R", myinp))
    }
  }
  return(shcommand) 
}

