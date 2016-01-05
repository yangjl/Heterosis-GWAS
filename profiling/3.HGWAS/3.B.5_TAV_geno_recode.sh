

#####
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





