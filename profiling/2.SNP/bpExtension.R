### Jinliang Yang
### Updated: 6/5/2012
### using GBS data to increase the resolution at the breakpoint of NAM RILs



bpExtension <- function(bpfile="bp_tagging1000.txt", GBSfile="tsnp_GBS50k_AGPv2.txt", parellel=TRUE, 
                        range=1:500, temfile="temfile0", outputfile="bp_GBS50k.txt"){
#########################################################################################################
# bpextension is a function to calculate the fine breakpoint of a given RIL by GBS snp data
# bpfile: breakpoint file
# GBSfile: GBS SNP from buckler lab
# parellel=TRUE, do the cal simulateously
# range: the range of cal each time
# temfile: the name of the temp file during the calculation.
  
  
  ### read in the breakpoint data ####
  bpfile <- read.table(bpfile, header=F)
  names(bpfile) <- c("id", "chr", "pos1","pos2", "snp1", "snp2");
  
  ### read in the first 5 rows of GBS data ###
  gbs5rows <- read.table(GBSfile, header = TRUE, nrows = 5)
  gbsnm <- names(gbs5rows)
  
  allrils <- as.character(unique(bpfile$id));
  if(parellel){
    allrils <- allrils[range];
  }

  count =0;

for(rilname in allrils){
  idx <- which(gbsnm == rilname)
  mynamril = subset(bpfile, id==rilname);
  
  if(length(idx) != 0){
  command = paste("awk '{print $1, $2, $3, $4, $5, $", idx, "}' tsnp_GBS50k_AGPv2.txt > ", temfile, sep="")
  system(command);
  gbs <- read.table(temfile, header=TRUE)
  
  
  
  for(i in 1:nrow(mynamril)){
    snp1 <- mynamril$snp1[i];
    snp2 <- mynamril$snp2[i];
    if(snp1 != snp2){
      subgbs <- subset(gbs, ch==mynamril$chr[i] & pos >= mynamril$pos1[i] & pos<mynamril$pos2[i] )
      subgbs0 <- subgbs[subgbs[,6]==0,];
      subgbs2 <- subgbs[subgbs[,6]==2,];
      
      if(nrow(subgbs2)>0){
        startmrk2 <- subgbs2[1,]
        idx1 <- which(subgbs$marker == startmrk2$marker)
        endmrk2 <- subgbs2[nrow(subgbs2),]
        idx2 <- which(subgbs$marker == endmrk2$marker)
      }
      ### initial the pos1 and pos2
      pos1 = mynamril$pos1[i];
      pos2 = mynamril$pos2[i];
      
      ### For cases of 0-----2 ###
      if(snp1 == 0){
        if(nrow(subgbs2) == 0 & nrow(subgbs0) > 0) {
          pos1 <- subgbs[nrow(subgbs0),]$pos;
          pos2 <- mynamril$pos2[i];
        }
        if(nrow(subgbs2) >0 & idx1 == 1){
          pos1 <- mynamril$pos1[i];
          pos2 <- startmrk2$pos;
        }
        if(nrow(subgbs2) >0 & idx1 > 1){
          pos1 <- subgbs[idx1-1,]$pos;
          pos2 <- subgbs[idx1,]$pos;
        } 
      }
      
      ### For cases of 2-----0 ###
      if(snp2 ==0){
        if(nrow(subgbs2) ==0 & nrow(subgbs0) >0){
          pos1 <- mynamril$pos1[i];
          pos2 <- subgbs0[1,]$pos;
        }
        if(nrow(subgbs2) >0 & idx2 == nrow(subgbs)){
          pos1 <- subgbs[idx2,]$pos;
          pos2 <- mynamril$pos2[i];
        }
        if(nrow(subgbs2) >0 & idx2 < nrow(subgbs)){
          pos1 <- subgbs[idx2,]$pos;
          pos2 <- subgbs[idx2+1,]$pos;
        }
      }
      
    #### assign pos1 and pos2 to the right pos ####
    mynamril$pos1[i] <- pos1;
    mynamril$pos2[i] <- pos2;
    mynamril$pos2[i-1] <- pos1;
    mynamril$pos1[i+1] <- pos2;
      if(pos2 < pos1){
        stop(paste(rilname, "at chr", mynamril$chr[i], pos1, pos2, "error!", sep=" "));
        
      }
    } 
  }
  }
  count = count+1;
  print(paste("###---Finished the", count, "RIL:", rilname,"---###", sep=" "));
  mynamril$dis <- mynamril$pos2-mynamril$pos1;
  
  write.table(mynamril, outputfile, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE )
}#End of the for loop ###
} #end of the function ####

