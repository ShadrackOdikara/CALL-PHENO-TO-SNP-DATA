library(dplyr)
library(nortest)
#LOAD THE DATASETS INTO R ENVIRONMENT

pheno <- read.table("PH_Analysed.txt", header = TRUE, sep = "\t")

hist(pheno$GeoM)	#show a distribution that;s a bit skewed to the left

ad.test(pheno$GeoM)	#Insufficient evidence to reject H0 of normality.

#Quick and dirty quatile normalization of the ear diameter data. Then, show
# histogram, and test for normality

GeoM_norm <- data.frame(Taxa = pheno$Taxa,
                      GeoM_norm = qqnorm
                      (pheno$GeoM, plot = F)$x)


hist(GeoM_norm$GeoM_norm)
ad.test(GeoM_norm$GeoM_norm)
head(GeoM_norm)
pheno1 <- cbind(GeoM_norm, pheno[,3])
head(pheno1)
#SELECT SPECIFIC COLUMNS IN THE SNP DATASET THAT MATCH THE AVALIABLE PHENOTYPE
library(data.table)
dt1 <- fread("ADP_genome.hmp.txt", 
            #select=c(pheno[,1]), 
            header = TRUE)
View(dt1)
dt <- fread("ADP_genome.hmp.txt", 
            select=c(pheno[,1]), 
            header = TRUE)
#View(dt)
#head(dt)
#dim(dt)
dt2 <- dt1[,1:11]
#dim(dt2)
#dt2
dim(dt)
geno <- cbind(dt2, dt)
View(geno)
dim(pheno1)

#View(geno)
dim(geno)
#write.csv(geno, file = "genotypes.csv")
geno1 <- read.csv("genotypes.csv", header = FALSE)
#write.ftable(geno, file = "genotype.csv")

#write.csv(pheno, file = "pheno.csv")
#Load the reqiured libraries (according to GAPIT documentation)
suppressPackageStartupMessages({
  library(multtest)
  library(gplots)
  library(LDheatmap)
  library(genetics)
  library(ape)
  library(EMMREML)
  library(compiler)
  library(scatterplot3d)
})




#############################################################
#############################################################

#load the GAPIT source code and helper function.
#source("gapit_function.R")
#source("emma.R")

#devtools::install_github("jiabowang/GAPIT3")

library(GAPIT3)


myGAPIT <- GAPIT(
  Y=pheno[,c(1,3)],
  G=geno1,
  CV=pheno[,c(1,2,4,5)],
  PCA.total = 7,
  #LD = 0.25,
  #LD.range = 500000,
  LD.chromosome = 500000,
  SNP.FDR = 1,
  SNP.MAF = 0.02,
  #MAF.calculate=TRUE,
  model = c("FarmCPU"),
  #p.threshold = 0.005,
  cutOff = 0.1,
  #DPP = 0.005,
  #alpha = c(0.05),
  Major.allele.zero = T,
  Inter.Plot = TRUE
  #h2 = 0.46,
  #SNP.effect = "Dom",
  #SNP.permutation = TRUE
)




