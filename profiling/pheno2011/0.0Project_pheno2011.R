# Jinliang Yang
# updated: 7/22/2012
# main script for the pheno2011 raw data processing

# Create the library:
library('ProjectTemplate')
create.project('pheno2011')


#####---------------data munging = data pre-processing----------------------------########

############################################
#####            Diallel                 ###
############################################

### pre-procssing the raw data###
source("munge/01.data_diallele_raw_processing.r")
cache('diallel.post.raw') # 3886 row x 22 cols with decription for each data point

##### BLUE of the 225 Diallels ####
diallel <- diallel.post.raw
source("munge/02-A.pheno_diallele_BLUE.r")

##### BLUE of the 27 NAM parents ####
source("munge/02-B.pheno_panzea_parents.r")

##### using BLUE to calculate HPH and MPH ####
source("munge/02-C.pheno_diallel_heterosis.r")

##### to calculate GCA and SCA ####
source("munge/02-D.pheno_diallel_SCA_SAS.r")

############################################
#####            2011 BxNAMRIL           ###
############################################
# For every script, get the long format and BLUE;

source("munge/03-A.data_NAMF1_KRN.r");
source("munge/03-B.data_NAMF1_cob.r");
source("munge/03-C.data_NAMF1_kernel.r");

#data merging
source("munge/03-D.data_NAMF1_merging.r")
cache('nam2011')
cache('nam2011blue');

############################################
#####            2008 IBM                ###
############################################
# For every script, get the long format and BLUE;
source("munge/04-A.IBM_data_processing.r")
cache('pheno2008')


############################################
#####         NAM from patrick brown     ###
############################################
# For every script, get the long format and BLUE;

source("04-B.NAM_pb_data1.r")
source("04-C.NAM_buckler_data.r")

cl <- rbind(pbcl, nam2011cl, cl2008);
cl <- cl[!is.na(cl$CL),]
cl_blue <- BLUE(data=cl, model=CL~Genotype, random=~1|Location, trait="CL", intercept="Z001E0001")

cd <- rbind(pbcd, nam2011cd, cd2008);
cd <- cd[!is.na(cd$CD),]
cd_blue <- BLUE(data=cd, model=CD~Genotype, random=~1|Location, trait="CD", intercept="Z001E0001")

nam_blue[['cl']] <- cl_blue;
nam_blue[['cd']] <- cd_blue;
cache('nam_blue')


#############################################
#####         Bx NAM RIL and IBM RILs     ###
#####         Mx NAM RILs                 ###
#############################################

source("04-D.NAM_Bx_Mx.r")
cache('BxMx') # is a list with two data.frames of Bx and Mx data





### test for consistency
inflo07 <- read.csv("~/Documents/VariationDB/pheno/INFLO7_BLUPs.csv")
test <- merge(cl_blue, inflo07, by.x="Genotype", by.y="line")
cor(test$CL.x, test$CL.y)
plot(test$CL.x, test$CL.y)















