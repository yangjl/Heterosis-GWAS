#Jinliang yang
#Purpose: Use python code to merge snps from 3 datasets
#updated: 7.2.2014

#Note: Running on server 129.186.85.7
run_snp3merge <- function(chr="chr1"){
  hmp1 <- paste("/mnt/02/yangjl/DBcenter/VariationDB/HapMap1/hapmap1V2_070214_", chr, ".dsf", sep="");
  hmp2 <- paste("/mnt/02/yangjl/DBcenter/VariationDB/HapMap2/maizeHapMap2V2_", chr, ".dsf", sep="");
  rnaseq <- paste("/mnt/02/yangjl/DBcenter/VariationDB/RNA-seq/rnaseq.AGPv2.condensed_", chr, ".dsf", sep="")
  file=data.frame(type=c("hapmap1", "hapmap2", "rnaseq"), dir=c(hmp1, hmp2, rnaseq));
  outfile <- paste("file_", chr, ".txt", sep="")
  write.table(file, outfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  out <- paste("snp_merged_", chr, sep="")
  command <- paste("python ~/Documents/PyCodes/snp3merge_v2.py -f", outfile, "-m 0", "-o", out);
  try(system(command))
}


setwd("~/DBcenter/VariationDB/merged/")
run_snp3merge(chr="chr1")
run_snp3merge(chr="chr2")
run_snp3merge(chr="chr3")
run_snp3merge(chr="chr4")
run_snp3merge(chr="chr5")
run_snp3merge(chr="chr6")
run_snp3merge(chr="chr7")
run_snp3merge(chr="chr8")
run_snp3merge(chr="chr9")
run_snp3merge(chr="chr10")


################ manually checked the merging results!
hapmap1 /mnt/02/yangjl/DBcenter/VariationDB/HapMap1/hapmap1V2_070214_chr1.dsf
hapmap2 /mnt/02/yangjl/DBcenter/VariationDB/HapMap2/maizeHapMap2V2_chr1.dsf
rnaseq /mnt/02/yangjl/DBcenter/VariationDB/RNA-seq/rnaseq.AGPv2.condensed_chr1.dsf

python ~/Documents/PyCodes/snp3merge_v2.py -f file_dir.txt -m 0 -o chr1_merged


##test the codes

setwd("~/DBcenter/VariationDB/merged/")
rnaseq <- read.table("~/DBcenter/VariationDB/RNA-seq/rnaseq.AGPv2.condensed_chr1.dsf", header=T)
hmp1 <- read.table("~/DBcenter/VariationDB/HapMap1/hapmap1V2_070214_chr1.dsf", header=T)
hmp2 <- read.table("~/DBcenter/VariationDB/HapMap2/maizeHapMap2V2_chr1.dsf", header=T)

subset(hmp1, snpid==snpid)
subset(hmp2, snpid==snpid)
subset(rnaseq, snpid==snpid)



