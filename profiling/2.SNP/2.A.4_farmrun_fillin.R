### Jinliang Yang
### use impute_parent in CJ data

#$SLURM_ARRAY_TASK_ID $SLURM_JOB_ID
source("~/Documents/Github/zmSNPtools/Rcodes/setUpslurm.R")

mysh <- paste("module load gcc jdk/1.8 tassel/5", 
"run_pipeline.pl -Xmx64g -fork1 -FILLINFindHaplotypesPlugin -hmp /home/jolyang/dbcenter/NAMP/namp_14m_maf5miss1.hmp.txt -o /home/jolyang/dbcenter/NAMP/haps -runfork1", sep="\n")

setUpslurm(slurmsh = "largedata/scripts/run_fillin.sh", 
           codesh = mysh, jobid = "run_fillin", wd=NULL, email="yangjl0930@gmail.com")


###>>> In this path: cd /home/jolyang/Documents/Github/Heterosis-GWAS
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --ntasks=1 largedata/scripts/run_fillin.sh

#$SLURM_ARRAY_TASK_ID $SLURM_JOB_ID
source("~/Documents/Github/zmSNPtools/Rcodes/setUpslurm.R")
mysh <- paste("module load gcc jdk/1.8 tassel/5", 
              "run_pipeline.pl -Xmx64g -fork1 -FILLINImputationPlugin -hmp /home/jolyang/dbcenter/AllZeaGBS/ZeaGBSv27_NAM_5873.hmp.txt -d /home/jolyang/dbcenter/NAMP/haps -o /home/jolyang/dbcenter/AllZeaGBS/ZeaGBSv27_NAM_5873_imputed.hmp.txt.gz -runfork1", sep="\n")

setUpslurm(slurmsh = "largedata/scripts/run_fillin.sh", 
           codesh = mysh, jobid = "run_fillin", wd=NULL, email="yangjl0930@gmail.com")

#$SLURM_ARRAY_TASK_ID $SLURM_JOB_ID
source("~/Documents/Github/zmSNPtools/Rcodes/setUpslurm.R")
mysh <- paste("module load gcc jdk/1.8 tassel/5", 
              "run_pipeline.pl -Xmx64g -fork1 -FILLINImputationPlugin -hmp /home/jolyang/dbcenter/AllZeaGBS/ZeaGBSv27_NAM_5873.hmp.txt -d /home/jolyang/dbcenter/NAMP/haps -o /home/jolyang/dbcenter/AllZeaGBS/ZeaGBSv27_NAM_5873_imputed.hmp.txt.gz -accuracy true -runfork1", sep="\n")

setUpslurm(slurmsh = "largedata/scripts/run_fillin2.sh", 
           codesh = mysh, jobid = "run_fillin2", wd=NULL, email="yangjl0930@gmail.com")
