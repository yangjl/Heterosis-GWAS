# Formatting the 2011 Diallel phenotype data
# Jinliang Yang
# Jan 17th, 2012
# updated: 4.6.2012
# server .7

##########################################################################################################################
# Check the data completeness and correct the notes
#Notes: 
#Dairy 8651-9000 v
#Zumwalt 9251-9600 v
#Johnson J2601-J2950 v
#Ki11: kill@ => Ki11 v
#Ms71: MS71 => Ms71 v

file1 <- read.csv("~/Documents/workingSpace/Heterosis_GWAS/pheno2011/2011rawdata/LL_KRN_test.csv")
dim(file1)
#8678 8

d1 <- subset(file1, Note.y=="Diallel" | Note.y=="NAM Filler" | Note.y=="B73");
dim(d1)
#[1] 3238    8
summary(d1)
############################################
### Delete duplicated input and missing data
# Get ride of 4 duplicated records
dup <- which(duplicated(d1$Barcode))
d1 <- d1[!duplicated(d1$Barcode),]

# Get ride of the missing records
### checking for outlyers
d1$KRN <- as.numeric(as.character(d1$KRN));
d1[is.na(d1$KRN),]
d1 <- d1[!is.na(d1$KRN),]
hist(d1$KRN, breaks=10)
dim(d1)
#[1] 3225    8

#setwd("/Users/yangjl/Documents/workingSpace/pheno2011/codes")
##########################################################################################################################
# Check the data completeness and correct the notes
#Notes: 
#Dairy 8651-9000 v
#Zumwalt 9251-9600 v
#Johnson J2601-J2950 v
#Ki11: kill@ => Ki11 v
#Ms71: MS71 => Ms71 v

file2 <- read.csv("data/LL_KN_AKW_test.csv")
d2 <- subset(file2, Note.y=="Diallel" | Note.y=="NAM Filler"| Note.y=="B73");
dim(d2)
#[1] 3204   11
summary(d2)

############################################
### Delete duplicated input and missing data
# Get ride of three duplicated records
barcode <- d2[duplicated(d2$Barcode),]$Barcode
d2[d2$Barcode %in% barcode,]
d2 <- d2[!duplicated(d2$Barcode),]

# Get ride of recodes with note.x nad KC<5
#d2 <- subset(d2, Note.x=="" & KC >= 50)
#d2[is.na(d2$TKW),]
#d2 <- d2[!is.na(d2$KRN),]
#########################################
### checking for outlyers
d2$KC <- as.numeric(as.character(d2$KC));
dim(d2)
# 3200 11
hist(d2$AKW)
hist(d2$KC, xlim=c(0, 1000))
hist(d2$TKW)

#setwd("/Users/yangjl/Documents/workingSpace/pheno2011/2011rawdata")
##########################################################################################################################
# Check the data completeness and correct the notes
#Notes: 
#Dairy 8651-9000 v
#Zumwalt 9251-9600 v
#Johnson J2601-J2950 v
#Ki11: kill@ => Ki11 v
#Ms71: MS71 => Ms71 v

file3 <- read.csv("data/LL_cob.test.csv")
dim(file3)
#[1] 8774   10| kernel [1] 8286   11 |KRN #8678 8

d3 <- subset(file3, Note.y=="Diallel" | Note.y=="NAM Filler" | Note.y=="B73");
dim(d3)
#[1] 2966   10 | [1] 2819   11 | [1] 2850    8
summary(d3)

############################################
### Delete duplicated input and missing data
# Get ride of 126 duplicated records
barcode <- d3[duplicated(d3$Barcode),]$Barcode
d3[d3$Barcode %in% barcode,]
d3 <- d3[!duplicated(d3$Barcode),]

# Get ride of the missing records
#diallel[is.na(diallel$TKW),]
#diallel <- diallel[!is.na(diallel$KRN),]
#########################################
### checking for outlyers
d3$CL <- as.numeric(as.character(d3$CL));
d3$CD <- as.numeric(as.character(d3$CD));
d3$CW <- as.numeric(as.character(d3$CW));
idx <- which.max(d3$CW)
d3[idx,]
d3[idx, ]$CW <- 28.99

hist(d3$CL);
hist(d3$CD);
hist(d3$CW);
dim(d3)
#3231 10

id <- rbind(d1[,c("Barcode", "Row", "Pedigree", "Genotype", "Farm", "Note.y")], d2[, c("Barcode", "Row", "Pedigree", "Genotype", "Farm", "Note.y")], 
            d3[, c("Barcode", "Row", "Pedigree", "Genotype", "Farm", "Note.y")])
id <- id[!duplicated(id$Barcode),]
names(id)[6] <- "Type"
#B73    Diallel NAM Filler 
#390       2452        398
dim(id)
# 3240 6

###########################################################################################################################################################################
id$p1 <- NA;
id$p2 <- NA;
id$Genotype <- as.character(id$Genotype);
for(i in 1:nrow(id)){
  tem <- unlist(strsplit(id$Genotype[i], split=" x "));
  tem <- sort(tem)
  id$p1[i] <- toupper(tem[1]);
  id$p2[i] <- toupper(tem[2]);
}
id$p1 <- gsub("@", "", id$p1)

parents <- read.table("data/nam_parents", header=TRUE)
b73 <- data.frame(pop="B73", pedigree="B73", parent="B73")
parents <- rbind(b73, parents)
parents$parent <- toupper(parents$parent)
parents <- parents[, c(1,3)]

diallel <- merge(id, parents, by.x="p1", by.y="parent")
diallel <- merge(diallel, parents, by.x="p2", by.y="parent", all.x=TRUE)

#########################################################
diallel2 <- diallel[, c("Barcode", "Row", "Pedigree", "Genotype","Farm", "Type", "p1", "p2", "pop.x", "pop.y")]
names(diallel2) <- c("Barcode", "Row", "Pedigree", "Genotype","Farm", "Type", "p1", "p2", "pop1", "pop2")
diallel2$Genotype2 <- paste(diallel2$pop1, diallel2$pop2, sep="x")

diallel3 <- merge(diallel2, d1[, c("Barcode", "KRN", "Note.x")], by="Barcode", all.x=T)
names(diallel3)[13] <- "Note_KRN"
diallel3 <- merge(diallel3, d2[, c("Barcode", "KC", "TKW", "AKW", "Note.x", "Hand.count")], all.x=T)
names(diallel3)[17:18] <- c("Note_ear", "Hand_Count")
diallel3 <- merge(diallel3, d3[, c("Barcode", "CL", "CW", "CD", "Note.x")], all.x=T)
names(diallel3)[22] <- "Note_cob"

############################
# dim: 3886, 22: all the diallel records with pedigree and annotations
##
diallel.post.raw <- diallel3 


#write.table(diallel3, "reports/pheno_diallel_master_raw_040612.csv", sep=",", row.names=FALSE, quote=FALSE)










