# Jinliang Yang
# 7/24/2012
# 


#### Bx #################################################################################
getBx <- function(trait="KRN"){
	idx <- grep("B73 x ", nam2011blue$Genotype);
	nam2011Bx <- nam2011blue[idx,]
	nam2011Bx$Genotype <- gsub(" ", "", nam2011Bx$Genotype);
	B1 <- nam2011Bx[, c("Genotype", trait)]
	
	lowtrait <- tolower(trait);
	B2 <- pheno2008[[lowtrait]];
	B2 <- B2[, c("Genotype", "BxRIL")];
	names(B2)[2] <- trait;
	B2$Genotype <- paste("B73x", B2$Genotype);
	B2$Genotype <- gsub(" ", "", B2$Genotype);
	
	Bx <- rbind(B1, B2);
}

krn <- getBx(trait="KRN")
akw <- getBx(trait="AKW")
tkw <- getBx(trait="TKW")
kc <- getBx(trait="KC")
cd <- getBx(trait="CD")
cw <- getBx(trait="CW")
cl <- getBx(trait="CL")

bx <- merge(krn, akw, by="Genotype", all=TRUE)
bx <- merge(bx, tkw, by="Genotype", all=TRUE)
bx <- merge(bx, kc, by="Genotype", all=TRUE)
bx <- merge(bx, cd, by="Genotype", all=TRUE)
bx <- merge(bx, cw, by="Genotype", all=TRUE)
bx <- merge(bx, cl, by="Genotype", all=TRUE)

#### Mx #################################################################################
getMx <- function(trait="KRN"){
	lowtrait <- tolower(trait);
	B2 <- pheno2008[[lowtrait]];
	B2 <- B2[, c("Genotype", "MxRIL")];
	names(B2)[2] <- trait;
	B2$Genotype <- paste("Mo17x", B2$Genotype);
	B2$Genotype <- gsub(" ", "", B2$Genotype);
	
	return(B2)
}

krn <- getMx(trait="KRN")
akw <- getMx(trait="AKW")
tkw <- getMx(trait="TKW")
kc <- getMx(trait="KC")
cd <- getMx(trait="CD")
cw <- getMx(trait="CW")
cl <- getMx(trait="CL")

mx <- merge(krn, akw, by="Genotype", all=TRUE)
mx <- merge(mx, tkw, by="Genotype", all=TRUE)
mx <- merge(mx, kc, by="Genotype", all=TRUE)
mx <- merge(mx, cd, by="Genotype", all=TRUE)
mx <- merge(mx, cw, by="Genotype", all=TRUE)
mx <- merge(mx, cl, by="Genotype", all=TRUE)

BxMx <- list()
BxMx[['Bx']] <- bx;
BxMx[['Mx']] <- mx;





