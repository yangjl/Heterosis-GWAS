# Jinliang Yang
# 7/23/23012
# purpose: merging the data to a single file

#merge the long format
nam <- rbind(nam2011krn[, c("Row","Barcode","Pedigree","Genotype", "Note.y","Farm")],
		nam2011cob[, c("Row","Barcode","Pedigree","Genotype", "Note.y","Farm")],
		nam2011ear[, c("Row","Barcode","Pedigree","Genotype", "Note.y","Farm")])
nam <- nam[order(nam$Barcode),]
nam <- nam[!duplicated(nam$Barcode),]

nam <- merge(nam, nam2011krn[, c("Barcode", "KRN", "Note.x")], all=TRUE);
names(nam)[8] <- "Note_KRN";

nam <- merge(nam, nam2011cob[, c("Barcode", "CL", "CD", "CW", "Note.x")], all=TRUE);
names(nam)[12] <- "Note_Cob";

nam <- merge(nam, nam2011ear[, c("Barcode", "KC", "TKW", "AKW", "Note.x")], all=TRUE);
names(nam)[16] <- "Note_Ear";

nam2011 <- nam

#merge the BLUE files
genotype <- c(namkrn[,1], namcl[,1], namcd[,1], namcw[,1], namakw[,1], namtkw[,1], namkc[,1]);
genotype <- genotype[!duplicated(genotype)];

nam <- data.frame(Genotype=genotype, trait=0);

nam <- merge(nam, namkrn, by="Genotype", all=TRUE)
nam <- merge(nam, namcl, by="Genotype", all=TRUE)
nam <- merge(nam, namcd, by="Genotype", all=TRUE)
nam <- merge(nam, namcw, by="Genotype", all=TRUE)
nam <- merge(nam, namakw, by="Genotype", all=TRUE)
nam <- merge(nam, namtkw, by="Genotype", all=TRUE)
nam <- merge(nam, namkc, by="Genotype", all=TRUE)
nam <- nam[, -2]

nam2011blue <- nam








