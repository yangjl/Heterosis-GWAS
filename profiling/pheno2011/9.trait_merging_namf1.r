# purpose: merge and vis the NAM F1 pheno
# Jinliang Yang
# 1.21.2012
# how to generate the . problem?

setwd("~/Documents/Heterosis_GWAS/pheno2011");
krn <- read.csv("nam_F1_KRN.csv");
kernel <- read.csv("nam_F1_Kernel.csv")
cob <- read.csv("nam_F1_cob.csv")

trait <- merge(krn, cob, by="Genotype", all=TRUE);
trait <- merge(trait, kernel, by="Genotype", all=TRUE)


geno <- trait[duplicated(trait$Genotype),]$Genotype
length(geno)
trait[trait$Genotype%in%geno,]

trait <- trait[!duplicated(trait$Genotype),]
density(trait$KRN)

par(mfrow=c(2,4))
# KRN
krn <- density(trait[!is.na(trait$Note) & trait$Note=="RIL",]$KRN) # returns the density data 
bxkrn <- density(trait[!is.na(trait$Note) & trait$Note=="Bx",]$KRN) # returns the density data
plot(krn, col="red", xlim=c(10,22), ylim=c(0,0.35), main="KRN", xlab="row number") # plots the results
lines(bxkrn, col="blue")

# Kernel: KC, TKW and AKW
kc <- density(trait[!is.na(trait$Note) & trait$Note=="RIL",]$KC)
bxkc <- density(trait[!is.na(trait$Note) & trait$Note=="Bx",]$KC)
plot(kc, col="red", ylim=c(0, 0.005), xlim=c(0, 800), main="KC", xlab="count")
lines(bxkc, col="blue")

tkw <- density(trait[!is.na(trait$Note) & trait$Note=="RIL",]$TKW)
bxtkw <- density(trait[!is.na(trait$Note) & trait$Note=="Bx",]$TKW)
plot(tkw, col="red", xlim=c(0, 250), main="TKW", xlab="weight(g)")
lines(bxtkw, col="blue")

akw <- density(trait[!is.na(trait$Note) & trait$Note=="RIL",]$AKW)
bxakw <- density(trait[!is.na(trait$Note) & trait$Note=="Bx",]$AKW)
plot(akw, col="red", ylim=c(0, 15), main="AKW", xlab="weight(g)")
lines(bxakw, col="blue")

# cob: CL, CW and CD
cl <- density(trait[!is.na(trait$Note) & trait$Note=="RIL",]$CL)
bxcl <- density(trait[!is.na(trait$Note) & trait$Note=="Bx",]$CL)
plot(cl, col="red", xlim=c(50,220), ylim=c(0,0.03), main="CL", xlab="length(mm)" )
lines(bxcl, col="blue")

cw <- density(trait[!is.na(trait$Note) & trait$Note=="RIL",]$CW)
bxcw <- density(trait[!is.na(trait$Note) & trait$Note=="Bx",]$CW)
plot(cw, col="red", ylim=c(0,0.08), main="CW", xlab="weight(g)")
lines(bxcw, col="blue")

cd <- density(trait[!is.na(trait$Note) & trait$Note=="RIL",]$CD)
bxcd <- density(trait[!is.na(trait$Note) & trait$Note=="Bx",]$CD)
plot(cd, col="red", ylim=c(0,0.25), main="CD", xlab="diameter(mm)")
lines(bxcd, col="blue")




