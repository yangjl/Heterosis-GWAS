#Jinliang yang
#Purpose: formatting the HapMap2 V2 SNPs for imputation
#start: 6.6.2012
#updated: 7.1.2014

chr2dsf <- function(chrfile="maizeHapMapV2_B73RefGenV2_201203028_chr2.hmp.txt", output="maizeHapMap2_chr2_s1"){
  #### this function read in the maizeHapMap2 file by chr and then do the following:
  # 1. [remove the B73 != Ref SNP] do not do this anymore!
  # 2. remove the duplicated SNP
  # 3. [change the coding of -/+ to A/T];
  # 4. change the IUPAC coding to "N"
  # 5. Output to DSF (density snp format)
  
  begTime <- Sys.time();
  # read in the first 5 rows of the inputfile
  tab5rows <- read.delim(chrfile, header = TRUE, nrows = 5)
  classes <- sapply(tab5rows, class)
  chr <- read.delim(chrfile, header = TRUE, colClasses = classes);
  print("###---Finished reading the file ---###");
    
  # NAM files
  nam <- read.csv("~/DBcenter/VariationDB/NAM_populations.csv", header=TRUE)
  namp <- paste(toupper(nam$parent), ".MZ", sep="");
  
  chr <- chr[, c("rs.", "alleles", "chrom", "pos", "B73.MZ", namp)];
  names(chr) <- c("rs", "alleles", "chr", "pos", "B73", as.character(nam$zid));
  tot <- nrow(chr); 
  
  ### remove the duplicated rows ###
  names(chr)[1] <- "snpid"
  chr$snpid <- paste(chr$chr, chr$pos, sep="_")
  chr <- chr[!duplicated(chr$snpid), ]
  flt1 <- nrow(chr);
  
  
  ### remove the one that B73 is missing #######
  #chr <- subset(chr, B73 != "N" );
    
  ### change the -/+ to A/T  ###
  #alleles <- chr$alleles;
  #chr[chr == "-"] <- "A";
  #chr[chr == "+"] <- "T";
  #chr$alleles <- alleles;
    
  ### re-coding the "Y", "R" ... to "N" allele ###
  for(i in 5:31){
    chr[chr[,i]!="A" & chr[,i]!="T" & chr[,i] !="C" 
        & chr[,i] !="G" & chr[,i]!="N" & chr[,i]!="+" & chr[,i]!="-", i] <- "N"
  }
  
  ######################## OUTPUT for different purpose ##################################
  #for imputation use
  message(sprintf("total input [ %s ], duplicated [ %s ], remaining [ %s ]", tot, tot-flt1, nrow(chr)))
  message("Note IUPAC changed to N!")
  
  write.table(chr, output, sep="\t", quote=FALSE, row.names=FALSE)
  print(Sys.time() - begTime);
}

