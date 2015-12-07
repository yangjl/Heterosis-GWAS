#Jinliang yang
#Purpose: Use PLINK to trim the correlated SNPs
#start: 1.23.2012
#updated: 7.2.2012

#Note: Running on server 129.186.85.7

# Formatting Function:###################################################
# 1. [remove the B73 != Ref SNP] do not do this anymore!
# 2. remove the duplicated SNP
# 3. [change the coding of -/+ to A/T];
# 4. [change the IUPAC coding to "N"] no IUPAC for this input!
# 5. Output to DSF (density snp format)
######################################################################

#chr2ped <- function(chr=1){
	
hmp1_2dsf <- function(output="hapmap1V2_070214.dsf"){
  
  ### read in the data.frame
  tab5rows <- read.delim("maizeHapMapV1_B73RefGenV2_20110309_ALL.hmp.txt", header=T, nrow=5)
  classes <- sapply(tab5rows, class)
  hapmap1 <- read.delim("maizeHapMapV1_B73RefGenV2_20110309_ALL.hmp.txt", header=TRUE, colClasses = classes)
  tot <- nrow(hapmap1)
  
  hapmap1<- hapmap1[, c(1:4, 12:38)]
  names(hapmap1) <- c("rs","alleles","chr","pos","B73","B97","CML103",
                      "CML228","CML247","CML277","CML322","CML333","CML52","CML69","HP301",
                      "IL14H","KI11","KI3","KY21","M162W","M37W","MO17","MO18W","MS71","NC350","NC358",
                      "OH43","OH7B","P39","TX303","TZI8")
  
  names(hapmap1) <- c("snpid","alleles","chr","pos","B73","Z001","Z002",
                      "Z003","Z004","Z005","Z006","Z007","Z008","Z009","Z010",
                      "Z011","Z012","Z013","Z014","Z015","Z016","Z017","Z018","Z019","Z020","Z021",
                      "Z022","Z023","Z024","Z025","Z026")
  
  ### remove the duplicated rows ###
  hapmap1$snpid <- paste(hapmap1$chr, hapmap1$pos, sep="_");
  hapmap1 <- hapmap1[!duplicated(hapmap1$snpid), ]
  flt1 <- nrow(hapmap1);
  
  ######################## OUTPUT 
  message(sprintf("total input [ %s ], duplicated [ %s ], remaining [ %s ]", tot, tot-flt1, nrow(hapmap1)))
  for(i in 1:10){
    chr <- subset(hapmap1, chr==i);
    myoutput <- paste(output, "_chr", i, ".dsf", sep="")
    write.table(chr, myoutput, sep="\t", quote=FALSE, row.names=FALSE)
  }  
}

################################################################################
setwd("~/DBcenter/VariationDB/HapMap1")
hmp1_2dsf(output="hapmap1V2_070214")




