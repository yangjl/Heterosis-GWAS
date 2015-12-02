# Jinliang Yang
# purpose: tagging SNPs for Imputation
# Last update: 1.30.2011

setwd("/Users/yangjl/Documents/GWAS2_KRN/SNP/tSNP");
#install.packages("XML")
library(XML)

mrk <- read.csv("NAM_map_20080419.Eddy.csv")
mrk$AGPv2 <- NA;

url1 <- "http://www.panzea.org/db/searches/webform/marker_search?database=panzea&marker_name_operator=like_c&marker_name=";
url2 <- paste("&marker_type_operator=%3D&marker_type=all&gene_locus_operator=%3D&gene_locus=&chromosome_operator=%3D&chromosome=all&",
"from_position_operator=%3E%3D&from_position=&to_position_operator=%3C%3D&to_position=&order_by=asc&order_by=ama.chr&",
"order_by=asc&order_by=ama.chr_start&order_by=asc&order_by=cm.name&order_by=asc&order_by=ama.marker_type&output_format",
"=html&submit=Submit#results", sep="")


for (i in 1:nrow(mrk)){
#for (i in 2:10){
	snpnm <- as.character(mrk$marker[i]);
	snpnm <- unlist(strsplit(snpnm, "/"))[1]
	url <- paste(url1, snpnm, url2, sep="")
	
	doc <- htmlParse(url)
	tables <- getNodeSet(doc, "//table")
	table <- tables[[8]]
	d <- xmlValue(table[[3]]);
	mrk$AGPv2[i] <- d;
	print(i);
}


mrk$id <- NA;
mrk$chr <- -9;
mrk$v2 <- NA;
for (i in 1:nrow(mrk)){
	tem <- unlist(strsplit(mrk$AGPv2[i], "\n"));
	mrk$id[i] <- tem[1];
	mrk$chr[i] <- tem[5];
	mrk$v2[i] <- tem[6];
	
}

mrk2 <- mrk[, -6];
mrk2$v2 <- gsub(",", "",mrk2$v2);
mrk2[mrk2$chr=="\302\240",]$chr <- NA;
mrk2[mrk2$v2=="\302\240",]$v2 <- NA;

write.table(mrk2, "NAM_map_AGPv2_panzea_012712.csv", sep=",", row.names=FALSE, quote=FALSE);
########################################################################################################

# ---------------------------------
# Format the tag SNP data for imputation:
map <- read.csv("NAM_map_AGPv2_panzea_012712.csv")
############ remove markers mapped to different chr
#map <- subset(map, ch==chr |is.na(chr))
#dim(map)

tag <- read.table("markergenotypes_namibm_062608_no_pheno.txt")
ttag <- t(tag)

ttag <- ttag[c(-2,-3),]
write.table(ttag, "markergenotype_tranposed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

mytag <- read.delim("markergenotype_tranposed", header=TRUE)

#oldtag <- read.csv("NAM_genos_mapping_20080103.formated.csv", header=TRUE)
mytag <- cbind(map, mytag)
mytag1 <- subset(mytag, !is.na(v2) & ch==chr);
mytag1 <- mytag1[, -5:-7];
names(mytag1)[5] <- "pos";
mytag <- mytag1

dim(mytag)
#[1] 1055 4897 [1] 1007 4897
mytag[, 6:4897] [mytag[, 6:4897]!=0 & mytag[, 6:4897]!=2] <- -9

# 0=A=B73 2=B=nonB73
write.table(mytag, "tsnp_4impute_AGPv2.txt", sep="\t", row.names=FALSE, quote=FALSE)




