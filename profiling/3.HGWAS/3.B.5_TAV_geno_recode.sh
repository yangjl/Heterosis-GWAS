### Jinliang Yang

##### 758 SNPs
cd largedata/SNP/
qctool -g geno_13M_chr1.bgen -og chr1_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr2.bgen -og chr2_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr3.bgen -og chr3_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr4.bgen -og chr4_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr5.bgen -og chr5_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr6.bgen -og chr6_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr7.bgen -og chr7_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr8.bgen -og chr8_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr9.bgen -og chr9_tavsnp.gen -incl-rsids TAV_snpid.txt
qctool -g geno_13M_chr10.bgen -og chr10_tavsnp.gen -incl-rsids TAV_snpid.txt

#############################################
source("lib/genReCode.R")

chr <- data.frame()
for(i in 1:10){
    gen <- read.table(paste0("largedata/SNP/chr", i, "_tavsnp.gen"), header=FALSE)
    message(sprintf("###>>> recode [ chr %s ]", i))
    out <- genReCode(gen=gen, code="code0123", outfile=NULL)
    chr <- rbind(chr, out)
}

write.table(chr, "largedata/SNP/TAV_recoded.csv", sep=",", row.names=FALSE, quote=FALSE)



##### 5764 SNPs
cd largedata/SNP/
qctool -g geno_13M_chr1.bgen -og chr1_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr2.bgen -og chr2_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr3.bgen -og chr3_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr4.bgen -og chr4_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr5.bgen -og chr5_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr6.bgen -og chr6_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr7.bgen -og chr7_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr8.bgen -og chr8_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr9.bgen -og chr9_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt
qctool -g geno_13M_chr10.bgen -og chr10_tavsnp_5764.gen -incl-rsids TAV_snpid_5764.txt

#############################################
source("lib/genReCode.R")

#chr <- data.frame()
for(i in 1:10){
    gen <- read.table(paste0("largedata/SNP/chr", i, "_tavsnp_5764.gen"), header=FALSE)
    message(sprintf("###>>> recode [ chr %s ]", i))
    out <- genReCode(gen=gen, code="code0123", outfile=NULL)
    write.table(out, paste0("largedata/SNP/TAV_recoded_5764_chr", i, ".csv"), sep=",", row.names=FALSE, quote=FALSE)
}


out <- data.frame()
for(i in 1:10){
    tem <- fread(paste0("largedata/SNP/TAV_recoded_5764_chr", i, ".csv"))
    message(sprintf("###>>> loaded [ chr %s ]", i))
    out <- rbind(out, tem)
}

out <- as.data.frame(out)
write.table(out, "largedata/SNP/TAV_recoded_N5756.csv", sep=",", row.names=FALSE, quote=FALSE)

