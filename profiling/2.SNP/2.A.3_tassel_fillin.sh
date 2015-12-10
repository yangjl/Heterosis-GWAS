### Jinliang Yang

# copy the following codes to shell
srun.x11 -p bigmemh --mem 64000 --ntasks=8 --nodelist=bigmem4
module load gcc jdk/1.8 tassel/5

# to export - whatever the hdf5 file is - if you want plink -'Plink' or 'Hapmap' for hmp
run_pipeline.pl -Xmx64g -fork1 -h5 ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5 -includeTaxaInfile /home/jolyang/dbcenter/AllZeaGBS/Taxa_nam_5873.txt -export ZeaGBSv27_NAM_5873.hmp -exportType Hapmap -runfork1


# FILLINFindHaplotypesPlugin
# must be hmp.txt
run_pipeline.pl -Xmx64g -fork1 -FILLINFindHaplotypesPlugin -hmp /home/jolyang/dbcenter/NAMP/namp_14m_maf5miss1.hmp.txt -o /home/jolyang/dbcenter/NAMP/haps -runfork1
run_pipeline.pl -Xmx64g -fork1 -FILLINImputationPlugin -hmp /home/jolyang/dbcenter/AllZeaGBS/ZeaGBSv27_NAM_5873.hmp -d /home/jolyang/dbcenter/NAMP/haps -o <outFile.hmp.txt.gz> -runfork1




