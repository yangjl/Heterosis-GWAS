## Jinliang Yang
## 8/12/2014
## detect pleiotropic effect QTLs

sqtl <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/sep_qtl_table.csv")
sqtl <- subset(sqtl, pop != "BxRIL" & pop != "MxRIL")
dim(sqtl)
#[1] 549  12
table(sqtl$trait)
table(sqtl$pop)

################
pleioQTL <- function(sqtl=sqtl){
  myg <- 1
  sqtl$pleio <- 0
  
  out <- data.frame()
  for(pi in unique(sqtl$pop)){
    subp <- subset(sqtl, pop == pi)
    for(chi in 1:10){
      subchr <- subset(subp, chr==chi)
      subchr <- subchr[order(subchr$LeftCI_pos),]
      
      if(nrow(subchr) >= 2){
        subchr$pleio[1] <- myg
        for(j in 2:nrow(subchr)){
          if(subchr$RightCI_pos[j-1] > subchr$LeftCI_pos[j] & !(subchr$trait[j] %in% subchr$trait[1:(j-1)])){
            subchr$pleio[j] <- myg
          }else{
            myg <- myg+1
            subchr$pleio[j] <- myg
          }
        }
      }#end of if
      out <- rbind(out, subchr)
      ## after each chromosome, myg should increase!
      myg <- myg +1
    }
  }
  
  myt <- table(out$pleio)
  
  ptab <- data.frame(myt)
  ptab <- ptab[-1,]
  ptab <- subset(ptab, Freq >1)
  ptab$pid <- 1:nrow(ptab)
  
  message(sprintf("QTL 2 traits [ %s ], 3 traits[ %s ], 4 traits [%s], >4 [%s]", 
                  nrow(subset(ptab, Freq==2)), nrow(subset(ptab, Freq==3)),
                  nrow(subset(ptab, Freq==4)), nrow(subset(ptab, Freq>4)) ))
  
  out2 <- merge(out, ptab, by.x="pleio", by.y="Var1")
  return(out2)
}

###
q1 <- pleioQTL(sqtl=sqtl)
#QTL 2 traits [ 60 ], 3 traits[ 19 ], 4 traits [7], >4 [1]


################
pleioQTL2 <- function(jqtl=jqtl){
  myg <- 1
  jqtl$pleio <- 0
  
  out <- data.frame()
  for(chi in 1:10){
    subchr <- subset(jqtl, chr==chi)
    subchr <- subchr[order(subchr$LeftCI_pos),]
    
    if(nrow(subchr) >= 2){
      subchr$pleio[1] <- myg
      for(j in 2:nrow(subchr)){
        if(subchr$RightCI_pos[j-1] > subchr$LeftCI_pos[j] & !(subchr$trait[j] %in% subchr$trait[1:(j-1)])){
          subchr$pleio[j] <- myg
        }else{
          myg <- myg+1
          subchr$pleio[j] <- myg
        }
      }
    }#end of if
    out <- rbind(out, subchr)
    ## after each chromosome, myg should increase!
    myg <- myg +1
  }
  
  myt <- table(out$pleio)
  
  ptab <- data.frame(myt)
  ptab <- ptab[-1,]
  ptab <- subset(ptab, Freq >1)
  ptab$pid <- 1:nrow(ptab)
  
  message(sprintf("QTL 2 traits [ %s ], 3 traits[ %s ], 4 traits [%s], >4 [%s]", 
                  nrow(subset(ptab, Freq==2)), nrow(subset(ptab, Freq==3)),
                  nrow(subset(ptab, Freq==4)), nrow(subset(ptab, Freq>4)) ))
  
  out2 <- merge(out, ptab, by.x="pleio", by.y="Var1")
  return(out2)
}

jqtl <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/joint_qtl_table.csv")
q2 <- pleioQTL2(jqtl=jqtl)
#QTL 2 traits [ 13 ], 3 traits[ 9 ], 4 traits [1], >4 [0]
write.table(q2, "~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S.table4_pleioQTL.csv",
            sep=",", row.names=FALSE, quote=FALSE)

