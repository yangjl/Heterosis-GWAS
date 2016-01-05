# 2008 IBM Phenotype for Heterotic traits
# Jinliang YANG
# 7/24/2012


###############################################################################
# Step 1: read in the raw data of Ames
###############################################################################
### IBM KRN data
leaf <- read.csv("~/Documents/GWAS2_KRN/pheno/IBM_leaf_traits_20101025_KSU.csv")
leaf$MO.Number <- gsub("O", "0", leaf$MO.Number)
leaf <- leaf[, 1:2]
names(leaf) <- c("Genotype", "ID")


load("data/rawdata/pheno2008.RData")

for(i in 1:7){
	pheno2008[[i]] <- merge(leaf, pheno2008[[i]], by="ID", all=T);
	pheno2008[[i]]$Location <- "Ames08";
	
	tem <- pheno2008[[i]]
	tem1 <- tem[is.na(tem$Genotype),]
	tem2 <- tem[!is.na(tem$Genotype),];
	tem1$Genotype <- paste("Z017E", tem1$ID, sep="");
	temjoin <- rbind(tem1, tem2);
	
	pheno2008[[i]] <- temjoin
}




