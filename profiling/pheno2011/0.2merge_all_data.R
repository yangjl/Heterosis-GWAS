# Jinliang Yang
# updated: 8/9/2012
# merging all data and output for GenSel use

# Create the library:
setwd("/Users/yangjl/Documents/Heterosis_GWAS/pheno2011")
library('ProjectTemplate')
load.project()

#### NAM RILs BULE from IBM RILS, pb data and nam2011########
summary(nam_blue); #list, including the IBM RIL

#### BxNAM RILs pheno2008 and nam2011########
summary(pheno2008); #list
head(nam2011blue)

#### MxNAM RILs pheno2008 ########
summary(pheno2008); #list

#### Diallel 2011 ########
head(pheno.diallel.master.BLUE)


combine_data <- function(trait="akw", ylim=c(0, 20), xlim=c(0,0), printfile=FALSE){

	TRAIT <- toupper(trait);	
	namril <- nam_blue[[trait]]
	namril$Fix <- 1;
	namril$pop <- "NAMRIL";
	namril$FID <- gsub("E.*$", "", namril$Genotype);
	namril <- namril[!is.na(namril[,2]),]
	
	bxril <- nam2011blue[, c("Genotype", TRAIT)];
	bxril$Genotype <- gsub("B73 x ", "B73x", bxril$Genotype);
	bxril <- bxril[grep("B73x", bxril$Genotype),]
	bxril$Fix <- 1;
	bxril$pop <- "BxRIL";
	bxril$FID <- gsub("E.*$", "", bxril$Genotype);
	bxril$FID <- gsub("B73x", "", bxril$FID);
	bxril <- bxril[!is.na(bxril[,2]),];
	
	mxril <- pheno2008[[trait]];
	mxril <- mxril[, c("Genotype", "MxRIL")]
	names(mxril)[2] <- TRAIT;
	mxril$Genotype <- paste("Mo17x", mxril$Genotype, sep="")
	mxril$Fix <- 1;
	mxril$pop <- "MxRIL";
	mxril$FID <- "Z017";
	mxril <- mxril[!is.na(mxril[,2]),];
	
	diallel <- pheno.diallel.master.BLUE[, c("Genotype", TRAIT)]
	diallel$Fix <- 1;
	diallel$pop <- "Diallel";
	diallel$FID <- "Z000";
	diallel <- diallel[!is.na(diallel[,2]),];
	
	output <- rbind(namril, bxril, mxril, diallel);
	cat(table(output$pop), "\n")
	
	#################################################
	if(!xlim[2]){
		xlim = range(output[,2])
	}
	plot(density(subset(output, pop=="NAMRIL")[,2] ), col="red", lwd=2, xlim=xlim, ylim=ylim, main="Density Plot", xlab=TRAIT)
	lines(density(subset(output, pop=="Diallel")[,2]), col="bisque4", lwd=2)
	lines(density(subset(output, pop=="BxRIL")[,2]), col="cadetblue4", lwd=2)
	lines(density(subset(output, pop=="MxRIL")[,2]), col="chocolate4", lwd=2)
	abline(v=pheno.namparents.BLUE[, TRAIT][1], col="blue", lty=2, lwd=2)
	temp <- legend("topright", legend =c(" ", " ", " ", " "), title="Population", text.width=strwidth("NAMRIL"),
               lty=1, lwd=3, col=c("red", "bisque4", "cadetblue4", "chocolate4"), xjust=1, yjust=1)
	text(temp$rect$left + temp$rect$w, temp$text$y,
       c("NAMRIL", "Diallel", "BxRIL", "MxRIL"), pos=2)

	##################################
	if(printfile){
		set.seed(1234);
		part = 0.8
		idx1 <- sample(1:nrow(namril), floor(nrow(namril)*part));
		idx2 <- sample(1:nrow(bxril), floor(nrow(bxril)*part));
		idx3 <- sample(1:nrow(mxril), floor(nrow(mxril)*part));
		idx4 <- sample(1:nrow(diallel), floor(nrow(diallel)*part));
		
		train_set <- rbind(namril[idx1,], bxril[idx2,], mxril[idx3,], diallel[idx4,])
		val_set <- rbind(namril[-idx1,], bxril[-idx2,], mxril[-idx3,], diallel[-idx4,])
		
		file1 <- paste(trait, "_GenSel_fullset.txt", sep="")
		file2 <- paste(trait, "_GenSel_training.txt", sep="")
		file3 <- paste(trait, "_GenSel_validation.txt", sep="")
		
		write.table(output, file1, sep="\t", row.names=FALSE, quote=FALSE)
		write.table(train_set, file2, sep="\t", row.names=FALSE, quote=FALSE);
		write.table(val_set, file3, sep="\t", row.names=FALSE, quote=FALSE)
		
		cat("Three files has been printed to the following directory:\n",
			getwd(), "\n",
			file1, "\n",
			file2, "\n",
			file3, "\n"
		)
	}
	return(output)
	
}

akw <- combine_data(trait="akw", ylim=c(0, 20), xlim=c(0,0), printfile=TRUE)
tkw <- combine_data(trait="tkw", ylim=c(0, 0.02), xlim=c(0,0), printfile=TRUE)
kc <- combine_data(trait="kc", ylim=c(0, 0.015), xlim=c(0,0), printfile=TRUE)

cw <- combine_data(trait="cw", ylim=c(0, 0.1), xlim=c(0,0), printfile=TRUE)
cd <- combine_data(trait="cd", ylim=c(0, 0.4), xlim=c(0,0), printfile=TRUE)
cl <- combine_data(trait="cl", ylim=c(0, 0.03), xlim=c(0,0), printfile=TRUE)




