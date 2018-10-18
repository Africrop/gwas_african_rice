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
# Written by Philippe Cubry, Nhung Ta and Yves Vigouroux
#
###################################################################################################################################

# This script file was developed to produce plots for the association analysis
# It will also produce tables of candidates loci using a FDR threshold
# except for MLMM where the Bonferroni threshold suggested by the authors is used
# Information about arguments :
# the variable to test

#### Getting the arguments from the script invocation ####
args <- commandArgs(TRUE)
pheno_variable <- args[1]
ifelse(is.na(args[2]),fdr <- 5,fdr <- as.numeric(args[2]))
print(paste("FDR is set to",fdr))

# Defining libraries paths
.libPaths(c("/usr/local/R-3.3.3/lib64/R/library", .libPaths()))

print("loading libraries")
library(RColorBrewer); cols <- brewer.pal(9,"Set1")
library(qqman)
library(Cairo)
library(qvalue)

#### Function to calculate correlations and adding them to a plot ####
panel.cor <- function(x, y, digits = 5, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.00001) txt2 <- paste("p= ", "<0.00001", sep = "")
  text(0.5, 0.4, txt2)
}


#### Reading data ####
print("reading data")
data <- read.table(paste(pheno_variable,"/results_",pheno_variable,".txt",sep=""))
mlmm.threshold <- read.table(paste(pheno_variable,"/mlmm_threshold_",pheno_variable,".txt",sep=""))

#### Computing FDR ####
print("computing FDR")
for(a in c("lfmm","cate","emma","gapit")){
    assign(paste("qv.",a,".",pheno_variable,sep=""),
           qvalue(p = data[paste("pv.",a,sep="")],
                  fdr = fdr/100))
    assign(paste("candidates.",a,".",pheno_variable,sep=""),
           which(get(paste("qv.",a,".",pheno_variable,sep=""))$significant))
  }

#### Manhattan plots ####

##### PDF plots #####
print("producing PDF manhattan plots")
pdf(file = paste(pheno_variable,"/manhattan_plots_",pheno_variable,"_by_method.pdf",sep=""),paper = 'a4r')
par(mfrow=c(1,2))
#LFMM
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.lfmm",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.lfmm.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.lfmm.",pheno_variable,sep="")),
                                                        paste("pv.lfmm",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"LFMM"))
qq(pvector = data$pv.lfmm,main=paste("qqplot LFMM",pheno_variable))
#CATE
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.cate",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.cate.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.cate.",pheno_variable,sep="")),
                                                        paste("pv.cate",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"cate"))
qq(pvector = data$pv.cate,main=paste("qqplot cate",pheno_variable))
#EMMA
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.emma",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.emma.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.emma.",pheno_variable,sep="")),
                                                        paste("pv.emma",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"emma"))
qq(pvector = data$pv.emma,main=paste("qqplot emma",pheno_variable))
#GAPIT
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.gapit",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.gapit.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.gapit.",pheno_variable,sep="")),
                                                        paste("pv.gapit",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"GAPIT/MLM"))
qq(pvector = data$pv.gapit,main=paste("qqplot GAPIT/MLM",pheno_variable))

par(mfrow=c(1,1))
pairs(data[,c("pv.lfmm","pv.cate","pv.emma","pv.gapit")], upper.panel = panel.cor)

dev.off()

#### PNG format #####
print("producing PNG manhattan plots")
CairoPNG(filename = paste(pheno_variable,"/manhattan_plots_",pheno_variable,".png",sep=""),units = "cm",width = 21,height = 29.7,res = 600)
par(mfrow=c(4,2))
#LFMM
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.lfmm",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.lfmm.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.lfmm.",pheno_variable,sep="")),
                                                        paste("pv.lfmm",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"LFMM"))
qq(pvector = data$pv.lfmm,main=paste("qqplot LFMM",pheno_variable))
#CATE
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.cate",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.cate.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.cate.",pheno_variable,sep="")),
                                                        paste("pv.cate",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"cate"))
qq(pvector = data$pv.cate,main=paste("qqplot cate",pheno_variable))
#EMMA
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.emma",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.emma.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.emma.",pheno_variable,sep="")),
                                                        paste("pv.emma",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"emma"))
qq(pvector = data$pv.emma,main=paste("qqplot emma",pheno_variable))
#GAPIT
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.gapit",snp = "SNP",
          genomewideline = ifelse(test = length(get(paste("candidates.gapit.",pheno_variable,sep="")))>0,
                                  yes = -log10(max(data[get(paste("candidates.gapit.",pheno_variable,sep="")),
                                                        paste("pv.gapit",sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"GAPIT/MLM"))
qq(pvector = data$pv.gapit,main=paste("qqplot GAPIT/MLM",pheno_variable))

dev.off()

#### Correlation plots between associations methods #####
print("producing correlations plots between association methods")
CairoPNG(filename = paste(pheno_variable,"/corr_plots_",pheno_variable,".png",sep=""),units = "cm",width = 21,height = 29.7,res = 600)
par(mfrow=c(1,1))
pairs(data[,c("pv.lfmm","pv.cate","pv.emma","pv.gapit")], upper.panel = panel.cor)
dev.off()


#### Diagnostic qqplots #####
print("producing qq-plots")
CairoPNG(filename = paste(pheno_variable,"/diagnostic_plots_",pheno_variable,".png",sep=""),units = "cm",width = 21,height = 29.7,res = 600)
par(mfrow=c(4,2))
#AOV
plot(y = sort(runif(n = length(data$pv.aov))),
    x = sort(data$pv.aov),
     ylab="Expected p-values",
     xlab="Observed p-values",
     main="qqplot for Analysis of Variance") +
  abline(0,1,col="red")
#LFMM
plot(y=sort(runif(n = length(data$pv.lfmm))),
     x=sort(data$pv.lfmm),
     ylab="Expected p-values",
     xlab="Observed p-values",
     main="qqplot for LFMM") +
  abline(0,1,col="red")
#CATE
plot(y=sort(runif(n = length(data$pv.cate))),
     x=sort(data$pv.cate),
     ylab="Expected p-values",
     xlab="Observed p-values",
     main="qqplot for CATE") +
  abline(0,1,col="red")
#EMMA
plot(y=sort(runif(n = length(data$pv.emma))),
     x=sort(data$pv.emma),
     ylab="Expected p-values",
     xlab="Observed p-values",
     main="qqplot for EMMA") +
  abline(0,1,col="red")
#GAPIT
plot(y=sort(runif(n = length(data$pv.gapit))),
     x=sort(data$pv.gapit),
     ylab="Expected p-values",
     xlab="Observed p-values",
     main="qqplot for GAPIT/MLM") +
  abline(0,1,col="red")
#MLMM extBIC
plot(y=sort(runif(n = length(data$pv.mlmm.extBIC))),
     x=sort(data$pv.mlmm.extBIC),
     ylab="Expected p-values",
     xlab="Observed p-values",
     main="qqplot for MLMM/extBIC") +
  abline(0,1,col="red")
#MLMM mbonf
plot(y=sort(runif(n = length(data$pv.mlmm.mbonf))),
     x=sort(data$pv.mlmm.mbonf),
     ylab="Expected p-values",
     xlab="Observed p-values",
     main="qqplot for MLMM/mbonf") +
  abline(0,1,col="red")

dev.off()

#### MLMM plots ####
print("producing MLMM plots")
CairoPNG(filename = paste(pheno_variable,"/mlmm_",pheno_variable,".png",sep=""),units = "cm",width = 21,height = 29.7,res = 600)
par(mfrow=c(2,2))

#For extBIC
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.mlmm.extBIC",snp = "SNP",
          genomewideline = mlmm.threshold$V1,
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"MLMM extBIC"))
qq(pvector = data$pv.mlmm.extBIC,main=paste("qqplot MLMM extBIC",pheno_variable))

#For mbonf
manhattan(x = data,chr = "Chromosome.gapit",bp = "Pos",p = "pv.mlmm.mbonf",snp = "SNP",
          genomewideline = mlmm.threshold$V1,
          suggestiveline = FALSE, col=cols,main=paste(pheno_variable,"MLMM mbonf"))
qq(pvector = data$pv.mlmm.mbonf,main=paste("qqplot MLMM mbonf",pheno_variable))

dev.off()

#### Saving candidates lists ####
print("saving candidates lists")
for(a in c("lfmm","cate","emma","gapit")){
    write.table(x=(data[get(paste("candidates.",a,".",pheno_variable,sep="")),
                        c("SNP","Chr","Pos",paste("pv.",a,sep=""))]),
                file = paste0(pheno_variable,"/candidates_",pheno_variable,"_",a,"_fdr",fdr,"perc.txt"),
                quote = FALSE,
                row.names = FALSE)
}
write.table(x=data[which(-log10(data$pv.mlmm.extBIC) > mlmm.threshold$V1),
                      c("SNP","Chr","Pos","pv.mlmm.extBIC")],
            file = paste0(pheno_variable,"/candidates_",pheno_variable,"_mlmm_extBIC.txt"),
            quote = FALSE,
            row.names = FALSE)
write.table(x=data[which(-log10(data$pv.mlmm.mbonf) > mlmm.threshold$V1),
                      c("SNP","Chr","Pos","pv.mlmm.mbonf")],
            file = paste0(pheno_variable,"/candidates_",pheno_variable,"_mlmm_mbonf.txt"),
            quote = FALSE,
            row.names = FALSE)
print("Done")
