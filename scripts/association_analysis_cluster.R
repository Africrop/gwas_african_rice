###################################################################################################################################
#
# Copyright 2017 IRD and Grenoble-Alpes University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD and Grenoble-Alpes University
#
# Written by Philippe Cubry, Helene Pidon, Olivier Fran?ois and Yves Vigouroux
#
###################################################################################################################################

# This script file was developed to perform association analysis in African rice on the IRD bioinformatic cluster
# Information about arguments :
# first one need to be the phenotypic file name,
# second the variable to test

# For phenotypic data, please follow the following format -> file with column names,
# first need to be "Code_ADN" and list individuals names,
# following can be the phenotypic data, with missing data coded as "NaN"

# WARNING : only take complete matrix of genotypes


# Getting the arguments from the script invocation
args <- commandArgs(TRUE)

# Assigning variables
pheno_file <- args[1]
pheno_variable <- args[2]

# Defining libraries paths
.libPaths(c("/usr/local/R-3.3.3/lib64/R/library", .libPaths()))

# Launching libraries and functions
#source("http://www.bioconductor.org/biocLite.R")
library(LEA)
library(lfmm)
library(cate)
library(data.table)
library("gplots")
library("LDheatmap")
library("genetics")
library("MASS")
library("compiler")
library("RColorBrewer")
library("multtest")
library("scatterplot3d")
source("emma.txt") # modification of emma library
source("gapit_functions.txt") # source code for GAPIT
library("mlmm")

# Getting genotypic data
print("read genotypic data")
genotype <- fread("glab_genotype_imputed_with_metadata.txt",header = TRUE)

# Getting phenotypic data
print("read phenotypic data")
pheno <- read.table(pheno_file, na.strings = "NaN", header = TRUE)
pheno_data <- na.exclude(pheno[,c(1,which(colnames(pheno)==pheno_variable))])
row.names(pheno_data) <- pheno_data$Code_ADN

# Filtering genotypes according to accessions present in the phenotypic data
genotype <- genotype[,c(1,2,which(colnames(genotype) %in% row.names(pheno_data))),with=FALSE]

# Control that genotype and phenotype data are in the same order
print("ensure genotypes and phenotypes are in the same order")
genotype <- cbind(genotype[,c(1,2)],setcolorder(genotype[,-c(1,2)],row.names(pheno_data)))

# Filtering using MAF
## MAF function
maf = function(x){ 
  x = x[x != 9]
  f = sum(x) ; (min(f , length(x) - f ))/length(x)}

## Computing MAF per SNP
print("Compute MAF per SNP")
maf_genotype <- apply(genotype[,3:ncol(genotype)],1,maf)
genotype$maf <- maf_genotype

## Filtering data for MAF > 5%
genotype <- genotype[which(genotype$maf>0.05),]

## Filtering data for counfonding factors estimation
genotype.test <- genotype[which(genotype$maf>0.2),-c(1,2)]
genotype.test$maf <- NULL
genotype.test <- as.matrix(genotype.test)

# Saving list of retained SNPs with their maf in the corresponding dataset and list of individuals
retained_snps <- data.frame(genotype[["#CHROM"]],genotype[["POS"]],genotype[["maf"]])
genotype$maf <- NULL #remove the maf column
write.table(retained_snps,file = paste(pheno_variable,"/chrpos_",pheno_variable,".txt",sep=""))
individuals=colnames(genotype)[3:(ncol(genotype))]
write.table(individuals,file = paste(pheno_variable,"/indiv_",pheno_variable,".txt",sep=""))

# Creating intermediate files for GAPIT and lfmm analyses
print("creating required files for analyses")
genotype.gapit <- cbind(paste(genotype$`#CHROM`,genotype$POS,sep="_"),genotype[,3:ncol(genotype)])
genotype <- genotype[,3:ncol(genotype)]
genotype <- t(genotype)
write.table(pheno_data[2],paste(pheno_variable,"/pheno_",pheno_variable,".env",sep=""),row.names = FALSE,col.names = FALSE)
write.table(pheno_data,paste(pheno_variable,"/pheno2_",pheno_variable,".env",sep=""),row.names = FALSE,col.names = TRUE)
colnames(genotype) <- paste(retained_snps$genotype....CHROM...,retained_snps$genotype...POS...,sep="_")

# Performing the different kind of association analysis
## Analysis of Variance (AoV)
print("Performing ANOVA")
pv.aov<-vector(length = ncol(genotype))
for (i in 1:ncol(genotype)) {
  x<-genotype[, i]
  x<-unlist(x)
  res<-aov(pheno_data[,2]~as.factor(x))
  pv.aov[i]<-summary(res)[[1]][["Pr(>F)"]][[1]]
} 

## lfmm
print("Performing LFMM")
x = read.table(paste(pheno_variable,"/pheno_",pheno_variable,".env",sep=""))[,1]
x <- as.matrix(x)
factors.lfmm = lfmm_ridge( Y = t(genotype.test), X =scale(x), K = 4)
mod.lfmm <- lfmm_test( Y = as.matrix(genotype), X =scale(x), lfmm = factors.lfmm, calibrate = "gif")
pv.lfmm = mod.lfmm$calibrated.pvalue

## cate
print("Performing cate")
x = read.table(paste(pheno_variable,"/pheno_",pheno_variable,".env",sep=""))[,1]
mod.cate = cate( ~ x, X.data = data.frame(x), Y = as.matrix(genotype), r = 4, fa.method = "pc", calibrate = TRUE)
pv.cate = mod.cate$beta.p.value

## EMMA
print("Performing EMMA")
Kinship <- emma.kinship(t(genotype),"additive","all")   #genotypes imputés
mod.emma <- emma.REML.t( x, t(genotype), Kinship)  #x = phenotypic variable
pv.emma <- mod.emma$ps

## GAPIT
print("Performing GAPIT")
snps <- genotype.gapit[[1]] ; genotype.gapit[[1]] <- NULL
genotype.gapit <- t(genotype.gapit) ; colnames(genotype.gapit) <- snps
genotype.gapit[genotype.gapit==1] <- 2
myGD=data.frame(row.names(genotype.gapit),genotype.gapit) ; row.names(myGD) <- NULL
retained_snps<-read.table(paste(pheno_variable,"/chrpos_",pheno_variable,".txt",sep=""))[,c(1,2)]
myGM=cbind(colnames(myGD[,-1]),retained_snps)
colnames(myGM)=c("name","CHROM","POS")
myGM$CHROM=gsub("Chr","",myGM$CHROM)
myP=read.table(paste(pheno_variable,"/pheno2_",pheno_variable,".env",sep=""),header=TRUE,na.string="NaN")
myGAPIT<-GAPIT(Y=myP,GD=myGD,GM=myGM,PCA.total=4,memo="MLM")

## MLMM
print("Performing MLMM")
# Computing relationship matrix
average <- colMeans(genotype, na.rm = TRUE)
stdev <- apply(genotype,2,sd)
genotype.stand <- sweep(sweep(genotype,2,average,"-"),2,stdev,"/")
K_mat <- (genotype.stand %*% t(genotype.stand))/ncol(genotype.stand)
#Perform mlmm analysis
mod.mlmm <- mlmm(Y = pheno_data[,2],X = genotype,K = K_mat,maxsteps = 20,nbchunks = 4)
names(myGM) <- c("SNP","Chr","Pos");myGM$Chr <- as.integer(myGM$Chr)

## Setting directory into the variable folder
print("Gather and save the results")
setwd(paste("./",pheno_variable,sep=""))
# Saving outcomes of the analysis
chrpos <- read.table(paste("chrpos_",pheno_variable,".txt",sep=""))
chrpos$SNP <- paste(chrpos$genotype....CHROM...,chrpos$genotype...POS...,sep="_")
##Gather results from multiple analyses
results <- cbind(chrpos,pv.aov,pv.lfmm,pv.cate,pv.emma) ; colnames(results) <- c("Chr","Pos","MAF","SNP","pv.aov","pv.lfmm","pv.cate","pv.emma")
##Adding GAPIT results
results <- merge(results,myGAPIT$GWAS[c("P.value","FDR_Adjusted_P-values","SNP","Chromosome","Position ","maf")],by.x = "SNP",by.y="SNP",sort = FALSE)
colnames(results) <- c("SNP","Chr","Pos","MAF","pv.aov","pv.lfmm","pv.cate","pv.emma","pv.gapit","pv.gapit.fdr_adjusted","Chromosome.gapit","Position.gapit ","maf.gapit")
##Adding MLMM results
results <- merge(results,mod.mlmm$opt_extBIC$out,by.x = "SNP",by.y="SNP",sort = FALSE)
results <- merge(results,mod.mlmm$opt_mbonf$out,by.x = "SNP",by.y="SNP",sort = FALSE)
colnames(results) <- c("SNP","Chr","Pos","MAF","pv.aov","pv.lfmm","pv.cate","pv.emma","pv.gapit","pv.gapit.fdr_adjusted","Chromosome.gapit","Position.gapit ","maf.gapit","pv.mlmm.extBIC","pv.mlmm.mbonf")
write.table(results,paste("results_",pheno_variable,".txt",sep=""))
write.table(x=mod.mlmm$bonf_thresh,paste("mlmm_threshold_",pheno_variable,".txt",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)
save(mod.mlmm,file = paste("object_mlmm_",pheno_variable,".Rdata",sep=""))