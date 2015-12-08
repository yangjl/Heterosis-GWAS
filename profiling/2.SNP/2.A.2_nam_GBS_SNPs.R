### Jinliang Yang
### 12/7/2015

allzea <- read.csv("data/AllZeaGBSv2.7_publicSamples_metadata20140411.csv")
dim(allzea)
#[1] 17280    18

table(allzea$Project)
table(allzea$GermplasmSet)

nam <- subset(allzea, GermplasmSet %in% "NAM")

nam$taxa <- paste(nam$DNASample, nam$LibraryPrepID, sep=":")

write.table(nam$taxa, "~/dbcenter/AllZeaGBS/Taxa_nam_5873.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# Note: space delimted => %s/\n/ /g in vim
