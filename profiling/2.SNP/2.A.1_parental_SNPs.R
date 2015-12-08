### Jinliang Yang
### re-visit merged parental SNPs

library("data.table", lib="~/bin/Rlib/")

### step1
getsnpinfo <- function(){
    out <- data.frame()
    for(i in 1:10){
        chr <- fread(paste0("~/dbcenter/NAMP/snp_merged_chr", i, ".info"), header=TRUE)
        chr <- as.data.frame(chr)
        out <- rbind(out, chr)
    }
    return(out)
}


### step2
getsnptb <- function(){
    out <- data.frame()
    for(i in 1:10){
        chr <- fread(paste0("~/dbcenter/NAMP/snp_merged_chr", i, ".dsf"), header=TRUE)
        chr <- as.data.frame(chr)
        out <- rbind(out, chr)
    }
    return(out)
}


### step3
dsf2hmp <- function(snpinfo, namsnp){
    
    info$alleles <- paste(info$major, info$minor, sep="/")
    info$strand <- "+"
    info$assembly <- "AGPv2"
    info$center <- "Panzea"
    info$panelLSID <- "nam"
    info$QCcode <- NA
    
    info <- info[, c("snpid", "alleles", "chr", "pos", "strand", "assembly", "center", "MAF", "missing", "panelLSID", "QCcode")]
    
    out <- merge(info, namsnp[, -2:-3], by="snpid")
    names(out)[1:11] <- c("rs", "alleles", "chrom", "pos", "strand", "assembly", "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
    out <- out[order(out$chrom, out$pos),]
    return(out)
}

###################################################
info <- getsnpinfo()
namsnp <- getsnptb()
out <- dsf2hmp(snpinfo, namsnp)

write.table(out, "~/dbcenter/NAMP/namp_14m_maf5miss1.hmp", sep="\t", row.names=FALSE, quote=FALSE)


