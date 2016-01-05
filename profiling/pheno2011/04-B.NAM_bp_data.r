# Jinliang Yang
# 7/24/2012
# 

###############################################################################
# NAM data1 from Patrick Brown, CL and CD
###############################################################################

#Note that we transformed BN and ERN.
#The raw data file is large. Here is a dropbox link:
#The abbreviations in this file are: db=TL, dc=SL, dg=BN,si=CL,sf=CD,so=ERN.

inflo7 <- read.table("~/Documents/VariationDB/pheno/maize_inflo7_rawdata.txt", header=TRUE)
pbcl <- subset(inflo7, trait=="si");
pbcd <- subset(inflo7, trait=="sf");

getgeno <- function(pbkrn=pbcl){
	pbkrn$Genotype <- paste("ZZ", 1000+pbkrn$pop, "EE", 10000+pbkrn$entry, sep="")
	pbkrn$Genotype <- gsub("Z1", "", pbkrn$Genotype)
	pbkrn$Genotype <- gsub("E1", "", pbkrn$Genotype)
	table(pbkrn$location);
	
	return(pbkrn)
	
}

####### CL ###################
pbcl <- getgeno(pbkrn=pbcl);
pbcl <- pbcl[, c("Genotype", "location", "value")]
names(pbcl) <- c("Genotype", "Location", "CL")


nam2011$Genotype <- gsub("@", "", nam2011$Genotype)
nam2011cl <- nam2011[, c("Genotype", "Farm", "CL")];
nam2011cl <- nam2011cl[!is.na(nam2011cl$CL),]
idx1 <- grep("B73", nam2011cl$Genotype)
nam2011cl <- nam2011cl[-idx1,]
names(nam2011cl) <- c("Genotype", "Location", "CL")

cl2008 <- pheno2008[['cl']]
cl2008 <- cl2008[, c("Genotype", "Location", "RIL")]
names(cl2008) <- c("Genotype", "Location", "CL")





####### CD ###################
pbcd <- getgeno(pbkrn=pbcd);
pbcd <- pbcd[, c("Genotype", "location", "value")]
names(pbcd) <- c("Genotype", "Location", "CD")


nam2011$Genotype <- gsub("@", "", nam2011$Genotype)
nam2011cd <- nam2011[, c("Genotype", "Farm", "CD")];
nam2011cd <- nam2011cd[!is.na(nam2011cd$CD),]
idx1 <- grep("B73", nam2011cd$Genotype)
nam2011cd <- nam2011cd[-idx1,]
names(nam2011cd) <- c("Genotype", "Location", "CD")

cd2008 <- pheno2008[['cd']]
cd2008 <- cd2008[, c("Genotype", "Location", "RIL")]
names(cd2008) <- c("Genotype", "Location", "CD")




