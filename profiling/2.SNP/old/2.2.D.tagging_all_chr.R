### Jinliang Yang
### Updated: 6/27/2012
### preparing SNPs for all the chrs

########## tagging SNPs for NAM RILs
tag1 <- read.table("tagging_gbs50k_chr1.txt", header=T)
tag2 <- read.table("tagging_gbs50k_chr2.txt", header=T)
tag3 <- read.table("tagging_gbs50k_chr3.txt", header=T)
tag4 <- read.table("tagging_gbs50k_chr4.txt", header=T)
tag5 <- read.table("tagging_gbs50k_chr5.txt", header=T)
tag6 <- read.table("tagging_gbs50k_chr6.txt", header=T)
tag7 <- read.table("tagging_gbs50k_chr7.txt", header=T)
tag8 <- read.table("tagging_gbs50k_chr8.txt", header=T)
tag9 <- read.table("tagging_gbs50k_chr9.txt", header=T)
tag10 <- read.table("tagging_gbs50k_chr10.txt", header=T)
tag <- rbind(tag1, tag2, tag3, tag4, tag5, tag6, tag7, tag8, tag9, tag10);
tag <- tag[order(tag$id, tag$chr, tag$pos1),]
write.table(tag, "tagging_gbs50k_all.txt", row.names=FALSE, quote=FALSE, sep="\t");

########## tagging SNPs for IBM RILs
itag1 <- read.table("IBM_impute/tagging_IBM136_chr1.txt", header=T)
itag2 <- read.table("IBM_impute/tagging_IBM136_chr2.txt", header=T)
itag3 <- read.table("IBM_impute/tagging_IBM136_chr3.txt", header=T)
itag4 <- read.table("IBM_impute/tagging_IBM136_chr4.txt", header=T)
itag5 <- read.table("IBM_impute/tagging_IBM136_chr5.txt", header=T)
itag6 <- read.table("IBM_impute/tagging_IBM136_chr6.txt", header=T)
itag7 <- read.table("IBM_impute/tagging_IBM136_chr7.txt", header=T)
itag8 <- read.table("IBM_impute/tagging_IBM136_chr8.txt", header=T)
itag9 <- read.table("IBM_impute/tagging_IBM136_chr9.txt", header=T)
itag10 <- read.table("IBM_impute/tagging_IBM136_chr10.txt", header=T)
itag <- rbind(itag1, itag2, itag3, itag4, itag5, itag6, itag7, itag8, itag9, itag10);
itag <- itag[order(itag$id, itag$chr, itag$pos1),]
write.table(itag, "IBM_impute/tagging_IBM136_all.txt", row.names=FALSE, quote=FALSE, sep="\t");

########## tagging SNPs for Bx RILs
btag1 <- read.table("BMxRILs/tagging_Bx692_chr1.txt", header=T)
btag2 <- read.table("BMxRILs/tagging_Bx692_chr2.txt", header=T)
btag3 <- read.table("BMxRILs/tagging_Bx692_chr3.txt", header=T)
btag4 <- read.table("BMxRILs/tagging_Bx692_chr4.txt", header=T)
btag5 <- read.table("BMxRILs/tagging_Bx692_chr5.txt", header=T)
btag6 <- read.table("BMxRILs/tagging_Bx692_chr6.txt", header=T)
btag7 <- read.table("BMxRILs/tagging_Bx692_chr7.txt", header=T)
btag8 <- read.table("BMxRILs/tagging_Bx692_chr8.txt", header=T)
btag9 <- read.table("BMxRILs/tagging_Bx692_chr9.txt", header=T)
btag10 <- read.table("BMxRILs/tagging_Bx692_chr10.txt", header=T)
btag <- rbind(btag1, btag2, btag3, btag4, btag5, btag6, btag7, btag8, btag9, btag10);
btag <- btag[order(btag$id, btag$chr, btag$pos1),]
write.table(btag, "BMxRILs/tagging_Bx692_all.txt", row.names=FALSE, quote=FALSE, sep="\t");

########## tagging SNPs for Mx RILs
mtag1 <- read.table("BMxRILs/tagging_Mx289_chr1.txt", header=T)
mtag2 <- read.table("BMxRILs/tagging_Mx289_chr2.txt", header=T)
mtag3 <- read.table("BMxRILs/tagging_Mx289_chr3.txt", header=T)
mtag4 <- read.table("BMxRILs/tagging_Mx289_chr4.txt", header=T)
mtag5 <- read.table("BMxRILs/tagging_Mx289_chr5.txt", header=T)
mtag6 <- read.table("BMxRILs/tagging_Mx289_chr6.txt", header=T)
mtag7 <- read.table("BMxRILs/tagging_Mx289_chr7.txt", header=T)
mtag8 <- read.table("BMxRILs/tagging_Mx289_chr8.txt", header=T)
mtag9 <- read.table("BMxRILs/tagging_Mx289_chr9.txt", header=T)
mtag10 <- read.table("BMxRILs/tagging_Mx289_chr10.txt", header=T)
mtag <- rbind(mtag1, mtag2, mtag3, mtag4, mtag5, mtag6, mtag7, mtag8, mtag9, mtag10);
mtag <- mtag[order(mtag$id, mtag$chr, mtag$pos1),]
write.table(mtag, "BMxRILs/tagging_Mx289_all.txt", row.names=FALSE, quote=FALSE, sep="\t");




