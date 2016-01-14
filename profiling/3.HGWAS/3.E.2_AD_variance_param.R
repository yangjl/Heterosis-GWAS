### Jinliang Yang
### Jan 5th, 2016
### get degree of dominance


gblup_allchr <- function(out_pwd, out_gpar="gparameter.dat", out_snpe="output_snpeff_ce.snpe",
                         genofile="largedata/SNP/geno_gblup.txt", phenofile="largedata/pheno/pheno_sample.txt",
                         mapfile = "largedata/SNP/geno_gblup_map.txt",
                         trait_col){ 
    
    geno <- read.table(genofile, header=TRUE)
    pheno <- read.table(phenofile, header=TRUE)
    
    out_gpar <- paste0(out_pwd, out_gpar)
    cat(
        "1000 #numer of iterations",
        "1.0 2.0 7.0 #starting values of Va, Vd and Ve",
        "1.0e-08 #tolerance level",
        paste(phenofile, "#phenotype file"),
        paste(trait_col, "#trait position in phenotype file"),
        "2 #number of fixed factors",
        "4 5 #positions of fixed factors in phenotype file",
        "0 #number of covariables",
        "0 #positions of covariables in phenotype file",
        paste(nrow(pheno), "#number of individuals genotyped"),
        "1 #number of chromosomes",
        
        paste(ncol(geno)-1, genofile, "#genotype file for each chromosome"),

        paste0(out_pwd, "output_greml_ce #output file for GREML"),
        paste0(out_pwd, "output_gblup_ce #output file for GBLUP"),
        paste("#def_Q", 1),
        paste("#use_ai_reml", 1),
        paste("#iter_ai_reml_start", 3),
        paste("#missing_phen_val", -999),
        paste("#map_file",  mapfile), #snpid chr position
        paste0("output_mrk_effect ", out_pwd, out_snpe),
        file=out_gpar, append=FALSE, sep="\n")
    message(sprintf("###>>> run this [ greml_ce %s > %s/%s ]", out_gpar, out_pwd, gsub("snpe", "log", out_snpe)))
}

##########
pheno <- read.table("largedata/pheno/pheno_sample.txt", header=TRUE)
#[1] "sampleid" "ID_2"     "missing"  "pop"      "FID"      "order"   
#[7] "CD"       "KRN"      "AKW"      "CL"       "CW"       "KC"      
#[13] "TKW"
for(i in 7:13){
    gblup_allchr(out_pwd="largedata/snpeff/",
              out_gpar= paste0("gp_", names(pheno)[i], ".dat"), 
              out_snpe= paste0(names(pheno)[i], "_snpeff_ce.snpe"),
              genofile="largedata/SNP/geno_gblup.txt", phenofile="largedata/pheno/pheno_sample.txt",
              mapfile = "largedata/SNP/geno_gblup_map.txt",
              trait_col=i)
}

###>>> run this [ greml_ce largedata/snpeff/gp_asi_perse.dat > largedata/snpeff//asi_perse_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_dtp_perse.dat > largedata/snpeff//dtp_perse_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_dts_perse.dat > largedata/snpeff//dts_perse_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_eht_perse.dat > largedata/snpeff//eht_perse_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_gy_perse.dat > largedata/snpeff//gy_perse_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_pht_perse.dat > largedata/snpeff//pht_perse_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_tw_perse.dat > largedata/snpeff//tw_perse_logff_ce.log ]

###>>> run this [ greml_ce largedata/snpeff/gp_asi_pBPHmin.dat > largedata/snpeff//asi_pBPHmin_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_dtp_pBPHmax.dat > largedata/snpeff//dtp_pBPHmax_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_dts_pBPHmax.dat > largedata/snpeff//dts_pBPHmax_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_eht_pBPHmax.dat > largedata/snpeff//eht_pBPHmax_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_gy_pBPHmax.dat > largedata/snpeff//gy_pBPHmax_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_pht_pBPHmax.dat > largedata/snpeff//pht_pBPHmax_logff_ce.log ]
###>>> run this [ greml_ce largedata/snpeff/gp_tw_pBPHmax.dat > largedata/snpeff//tw_pBPHmax_logff_ce.log ]