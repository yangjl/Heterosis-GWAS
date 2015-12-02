## Jinliang Yang
## 7/21/2014
## detect synergistic effect QTLs

sqtl <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/sep_qtl_table.csv")
sqtl <- subset(sqtl, pop != "BxRIL" & pop != "MxRIL")
dim(sqtl)

table(sqtl$trait)
table(sqtl$pop)

################
synQTL <- function(sqtl=sqtl, trait1="KRN", trait2="AKW"){
  sub1 <- subset(sqtl, trait==trait1)
  sub2 <- subset(sqtl, trait==trait2)
  
  out <- data.frame();
  syn <- 0;
  ### from QTL1 find the corresponding syn QTL
  for(i in 1:nrow(sub1)){
    qtl2 <- subset(sub2, chr==sub1$chr[i] & pop==sub1$pop[i])
    if(nrow(qtl2) > 0){
      p1 <- sub1$LeftCI_pos[i];
      p2 <- sub1$RightCI_pos[i];
      effect <- sub1$effect[i];
      for(j in 1:nrow(qtl2)){
        q1 <- qtl2$LeftCI_pos[j]
        q2 <- qtl2$RightCI_pos[j]
        if(effect*qtl2$effect[j] >0 & 
             ((p1 < q1 &  p2 > q1) | (p1 < q2 & p2 > q2) | (p1 > q1 & p2 < q2))){
          syn <- syn +1;
          out <- rbind(out, sub1[i,], qtl2[j, ])
        } 
      }
    }
  }
  message(sprintf("Syn QTLs [ %s ]", syn))
  return(out)
}

###
krn_akw <- synQTL(sqtl=sqtl, trait1="KRN", trait2="AKW")
#Syn QTLs [ 6 ]
krn_cl <- synQTL(sqtl=sqtl, trait1="KRN", trait2="CL")
#Syn QTLs [ 5 ]
cd_cl <- synQTL(sqtl=sqtl, trait1="CD", trait2="CL")
#Syn QTLs [ 8 ]

syn <- rbind(krn_akw, krn_cl, cd_cl)
write.table(syn, "~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S.table3.synqtl.csv", 
            sep=",", row.names=FALSE,  quote=FALSE)