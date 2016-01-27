### Jinliang Yang
### 01/26/2016



tav <- read.csv("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/StableN.tav758.csv")


fgs <- read.delim("~/DBcenter/AGPv2/ZmB73_5b_FGS_info.txt", header=TRUE)
fgs <- subset(fgs, is_canonical == "yes")

classical <- read.delim("~/DBcenter/AGPv2/ZmB73_5a_named_genes.txt", header=TRUE)
classical <- read.csv("data/classical_maize_genes_v2.csv")
names(classical)[9] <- "gene_id"


count_tav_in_genes <- function(fgs, tav, classical, binsize=10000){
    
    ### 
    fgs$bin <- paste(fgs$chromosome, round(fgs$transcript_start/binsize, 0), sep="_")
    fgs$bin <- gsub("chr", "", fgs$bin)
    tav$bin <- paste(tav$chr, round(tav$pos/binsize, 0), sep="_")
    sum(tav$bin %in% fgs$bin)
    
    ###
    cla <- subset(fgs, gene_id %in% classical$gene_id)
    #a1 <- sum(tav$bin %in% cla$bin)
    a1 <- subset(cla, bin %in% tav$bin)
    gid1 <- length(unique(a1$gene_id))
    
    noncla <- subset(fgs, !(gene_id %in% classical$gene_id))
    #a2 <- sum(tav$bin %in% noncla$bin)
    a2 <- subset(noncla, bin %in% tav$bin)
    gid2 <- length(unique(a2$gene_id))
    
    M <- as.table(rbind(c(gid1, nrow(cla) - gid1), c(gid2, nrow(noncla) - gid2)))
    print(chisq.test(M))
    dimnames(M) <- list(hits = c("tavhit", "non-tavhits"),
                        geneclass = c("Classical","Others"))
    print(M)
    
    outgene <- merge(tav, cla, by="bin", all.x=TRUE)
    outgene <- merge(outgene, classical, by="gene_id", all.x=TRUE)
    return(outgene)
}

res <- count_tav_in_genes(fgs, tav, classical, binsize=10000)
res <- count_tav_in_genes(fgs, tav, classical, binsize=20000)
res <- count_tav_in_genes(fgs, tav, classical, binsize=30000)
res <- count_tav_in_genes(fgs, tav, classical, binsize=40000)
res <- count_tav_in_genes(fgs, tav, classical, binsize=60000)
res <- count_tav_in_genes(fgs, tav, classical, binsize=70000)

res <- count_tav_in_genes(fgs, tav, classical, binsize=50000)
res <- count_tav_in_genes(fgs, tav, classical, binsize=100000)

#res <- count_tav_in_genes(fgs, tav, classical, binsize=1000000)


subset(res, !is.na(gene_id))

