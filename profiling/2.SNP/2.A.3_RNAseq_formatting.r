#Jinliang yang
#Purpose: Use PLINK to trim the correlated SNPs
#updated: 6.7.2012

# Formatting Function:###################################################
# Modified form the GWAS_Hapmap1_SNP_prunning.r
######################################################################

#####---------RUN the PERL program to reformat file
#---./SNPmatrix_condense.pl NAM.AGPv2.SNPs+INDELs.pooled-20120217.txt > NAM.AGPv2.SNPs+INDELs_condensed.dsf 



rnaseq_2dsf <- function(output="~/test.dsf"){
  tab5rows <- read.table("NAM.AGPv2.SNPs+INDELs.condensed", nrow=5, header=FALSE)
  classes <- sapply(tab5rows, class)
  classes[2] <- "character"
  RNAseq <- read.table("NAM.AGPv2.SNPs+INDELs.condensed", header = FALSE, colClasses = classes)
  
  names(RNAseq) <- c("snpid", "chr", "pos","B73","B97","CML103", "CML228","CML247", 
                     "CML277","CML322","CML333","CML52","CML69","HP301",
                     "IL14H","Ki11","Ki3","Ky21","M162W","M37W",
                     "Mo17","Mo18W","Ms71","NC350","NC358","Oh43",
                     "Oh7B","P39","Tx303","Tzi8")
  
  names(RNAseq) <- c("snpid", "chr","pos", "B73", "Z001","Z002", "Z003","Z004", 
                     "Z005","Z006","Z007","Z008","Z009","Z010",
                     "Z011","Z012","Z013","Z014","Z015","Z016",
                     "Z017","Z018","Z019","Z020","Z021","Z022",
                     "Z023","Z024","Z025","Z026")
  RNAseq$alleles <- RNAseq$B73
  RNAseq <- RNAseq[, c("snpid", "alleles", "chr","pos", "B73", "Z001","Z002", "Z003","Z004", 
                       "Z005","Z006","Z007","Z008","Z009","Z010",
                       "Z011","Z012","Z013","Z014","Z015","Z016",
                       "Z017","Z018","Z019","Z020","Z021","Z022",
                       "Z023","Z024","Z025","Z026")]
  #remove the missing chr ==Mt, Pt and UNKNOWN
  RNAseq <- subset(RNAseq, chr %in% 1:10);
  tot <- nrow(RNAseq)
  #[1] 4364987      31
  RNAseq <- RNAseq[!duplicated(RNAseq$snpid), ]
  flt1 <- nrow(RNAseq);
  
  RNAseq <- RNAseq[order(RNAseq$chr, RNAseq$pos),]
  
  ######################## OUTPUT 
  message(sprintf("total variants [ %s ], duplicated [ %s ], remaining [ %s ]", tot, tot-flt1, nrow(RNAseq)))
  for(i in 1:10){
    chr <- subset(RNAseq, chr==i);
    myoutput <- paste(output, "_chr", i, ".dsf", sep="")
    write.table(chr, myoutput, sep="\t", quote=FALSE, row.names=FALSE)
  }  
}


setwd("~/DBcenter/VariationDB/RNA-seq")
rnaseq_2dsf(output="rnaseq.AGPv2.condensed")








