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
# Written by Philippe Cubry and Yves Vigouroux
#
###################################################################################################################################

# Defining libraries paths
.libPaths(c("/usr/local/R-3.3.3/lib64/R/library", .libPaths()))

#Required libraries
library(RColorBrewer); cols <- brewer.pal(9,"Set1")
library(qqman)
library(qvalue)
library(Cairo)

#Defining functions
##Function to calculate correlations and adding them to a plot
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

##Function to combine p-values from multiple tests using Fisher conbination probabilities test
comb.fisher <- function(x=pvalues_to_combine){
  # Calculating the statistic
  stat.combined.value <- -2*(sum(sapply(x,log)))
  degree_of_freedom <- 2*length(x)
  pchisq(stat.combined.value,df=degree_of_freedom,lower.tail=FALSE)
}


##Function to combine p-values using Z-transform test
comb.stouffer <- function(x=pvalues_to_combine){
  # Calculating the statistic
  combined_z_score <- sum(qnorm(x))/sqrt(length(x))
  resulting_pval <- pnorm(combined_z_score)
  resulting_pval
}

##Function to combine p-values using weigthed Z-method
comb.weightedZ <- function(x=pvalues_to_combine){
  # Calculating the statistic
  combined_z_score <- sum((as.numeric(indiv_number)*qnorm(x)))/sqrt(sum(as.numeric(indiv_number)^2))
  resulting_pval <- pnorm(combined_z_score)
  resulting_pval
}

# Getting the arguments from the script invocation
args <- commandArgs(TRUE)
suppressWarnings(if(is.na(as.numeric(args[length(args)]))){
  pheno_variables <- args ; fdr <- 5
} else {pheno_variables <- args[-length(args)] ; fdr <- as.numeric(args[length(args)])})

print(paste("FDR is set to",fdr))

pval_list <- vector(mode = "list",length(pheno_variables))
indiv_number <- vector(mode = "character",length(pheno_variables))
names(pval_list) <- NA
names(indiv_number) <- NA
  ## making a loop to retrieve the different datasets and putting them into a list, also retrieve the number of individuals for each test
for (i in 1:length(pheno_variables)){
  pval_list[[i]] <- read.table(file = paste(pheno_variables[i],"/results_",pheno_variables[i],".txt",sep=""))
  names(pval_list)[i] <- pheno_variables[i]
  names(pval_list[[i]])[names(pval_list[[i]])=="pv.lfmm"] <- paste0("pv.lfmm",pheno_variables[i])
  names(pval_list[[i]])[names(pval_list[[i]])=="pv.cate"] <- paste0("pv.cate",pheno_variables[i])
  names(pval_list[[i]])[names(pval_list[[i]])=="pv.emma"] <- paste0("pv.emma",pheno_variables[i])
  names(pval_list[[i]])[names(pval_list[[i]])=="pv.gapit"] <- paste0("pv.gapit",pheno_variables[i])
  indiv_number[i] <- nrow(read.table(file = paste(pheno_variables[i],"/indiv_",pheno_variables[i],".txt",sep="")))
  names(indiv_number)[i] <- pheno_variables[i]
}

# Creating dataframes for each kind of analysis
pv.lfmm <- merge(pval_list[[1]][c("SNP","Chr","Pos","Chromosome.gapit",
                                  grep(pattern = "pv.lfmm", x = names(pval_list[[1]]),value=TRUE))],
                 pval_list[[2]][c("SNP",
                                  grep(pattern = "pv.lfmm", x = names(pval_list[[2]]),value=TRUE))],
                 by="SNP",sort = FALSE)
if(length(pheno_variables)>2){
  for(i in 3:length(pheno_variables)){
    pv.lfmm <- merge(pv.lfmm ,
                     pval_list[[i]][c("SNP",
                                      grep(pattern = "pv.lfmm", x = names(pval_list[[i]]),value=TRUE))],by="SNP",sort = FALSE)
  }
}

pv.cate <- merge(pval_list[[1]][c("SNP","Chr","Pos","Chromosome.gapit",
                                  grep(pattern = "pv.cate", x = names(pval_list[[1]]),value=TRUE))],
                 pval_list[[2]][c("SNP",
                                  grep(pattern = "pv.cate", x = names(pval_list[[2]]),value=TRUE))],
                 by="SNP",sort = FALSE)
if(length(pheno_variables)>2){
  for(i in 3:length(pheno_variables)){
    pv.cate <- merge(pv.cate ,
                     pval_list[[i]][c("SNP",
                                      grep(pattern = "pv.cate", x = names(pval_list[[i]]),value=TRUE))],by="SNP",sort = FALSE)
  }
}

pv.emma <- merge(pval_list[[1]][c("SNP","Chr","Pos","Chromosome.gapit",
                                  grep(pattern = "pv.emma", x = names(pval_list[[1]]),value=TRUE))],
                 pval_list[[2]][c("SNP",
                                  grep(pattern = "pv.emma", x = names(pval_list[[2]]),value=TRUE))],
                 by="SNP",sort = FALSE)
if(length(pheno_variables)>2){
  for(i in 3:length(pheno_variables)){
    pv.emma <- merge(pv.emma ,
                     pval_list[[i]][c("SNP",
                                      grep(pattern = "pv.emma", x = names(pval_list[[i]]),value=TRUE))],by="SNP",sort = FALSE)
  }
}

pv.gapit <- merge(pval_list[[1]][c("SNP","Chr","Pos","Chromosome.gapit",
                                   grep(pattern = paste0("pv.gapit",pheno_variables[1]), x = names(pval_list[[1]]),value=TRUE))],
                 pval_list[[2]][c("SNP",
                                  grep(pattern = paste0("pv.gapit",pheno_variables[2]), x = names(pval_list[[2]]),value=TRUE))],
                 by="SNP",sort = FALSE)
if(length(pheno_variables)>2){
  for(i in 3:length(pheno_variables)){
    pv.gapit <- merge(pv.gapit ,
                     pval_list[[i]][c("SNP",
                                      grep(pattern = paste0("pv.gapit",pheno_variables[i]), x = names(pval_list[[i]]),value=TRUE))],
                     by="SNP",sort = FALSE)
  }
}

# Calculating combined statistic
pv.lfmm$Combined_test_fisher <- apply(pv.lfmm[,grep(pattern = "pv.lfmm", x = names(pv.lfmm),value=TRUE)],1,comb.fisher)
pv.cate$Combined_test_fisher <- apply(pv.cate[,grep(pattern = "pv.cate", x = names(pv.cate),value=TRUE)],1,comb.fisher)
pv.emma$Combined_test_fisher <- apply(pv.emma[,grep(pattern = "pv.emma", x = names(pv.emma),value=TRUE)],1,comb.fisher)
pv.gapit$Combined_test_fisher <- apply(pv.gapit[,grep(pattern = "pv.gapit", x = names(pv.gapit),value=TRUE)],1,comb.fisher)

pv.lfmm$Combined_test_stouffer <- apply(pv.lfmm[,grep(pattern = "pv.lfmm", x = names(pv.lfmm),value=TRUE)],1,comb.stouffer)
pv.cate$Combined_test_stouffer <- apply(pv.cate[,grep(pattern = "pv.cate", x = names(pv.cate),value=TRUE)],1,comb.stouffer)
pv.emma$Combined_test_stouffer <- apply(pv.emma[,grep(pattern = "pv.emma", x = names(pv.emma),value=TRUE)],1,comb.stouffer)
pv.gapit$Combined_test_stouffer <- apply(pv.gapit[,grep(pattern = "pv.gapit", x = names(pv.gapit),value=TRUE)],1,comb.stouffer)

pv.lfmm$Combined_test_weightedZ <- apply(pv.lfmm[,grep(pattern = "pv.lfmm", x = names(pv.lfmm),value=TRUE)],1,comb.weightedZ)
pv.cate$Combined_test_weightedZ <- apply(pv.cate[,grep(pattern = "pv.cate", x = names(pv.cate),value=TRUE)],1,comb.weightedZ)
pv.emma$Combined_test_weightedZ <- apply(pv.emma[,grep(pattern = "pv.emma", x = names(pv.emma),value=TRUE)],1,comb.weightedZ)
pv.gapit$Combined_test_weightedZ <- apply(pv.gapit[,grep(pattern = "pv.gapit", x = names(pv.gapit),value=TRUE)],1,comb.weightedZ)



CairoPNG(filename = paste("Combined_tests/Correlograms_lfmm_",paste(pheno_variables,collapse = "_"),".png",sep=""),
    units = "cm",width = 21,height = 29.7,res = 600)
correlogram.lfmm <- pairs(pv.lfmm[,5:length(pv.lfmm)], upper.panel = panel.cor)
dev.off()

CairoPNG(filename = paste("Combined_tests/Correlograms_cate_",paste(pheno_variables,collapse = "_"),".png",sep=""),
    units = "cm",width = 21,height = 29.7,res = 600)
correlogram.cate <- pairs(pv.cate[,5:length(pv.cate)], upper.panel = panel.cor)
dev.off()

CairoPNG(filename = paste("Combined_tests/Correlograms_emma_",paste(pheno_variables,collapse = "_"),".png",sep=""),
    units = "cm",width = 21,height = 29.7,res = 600)
correlogram.emma <- pairs(pv.emma[,5:length(pv.emma)], upper.panel = panel.cor)
dev.off()

CairoPNG(filename = paste("Combined_tests/Correlograms_gapit_",paste(pheno_variables,collapse = "_"),".png",sep=""),
    units = "cm",width = 21,height = 29.7,res = 600)
correlogram.gapit <- pairs(pv.gapit[,5:length(pv.gapit)], upper.panel = panel.cor)
dev.off()

# Computing FDR
for(i in pheno_variables){
  for(a in c("lfmm","cate","emma","gapit")){
  assign(paste("qv.",a,".",i,sep=""), qvalue(p = get(paste("pv.",a,sep=""))[paste("pv.",a,i,sep="")], fdr = fdr/100))
  assign(paste("candidates.",a,".",i,sep=""), which(get(paste("qv.",a,".",i,sep=""))$significant))
}}

for(a in c("lfmm","cate","emma","gapit")){
  for(i in c("fisher","stouffer",'weightedZ')){
  assign(paste("qv.",a,".",i,sep=""), qvalue(p = get(paste("pv.",a,sep=""))[paste("Combined_test_",i,sep="")], fdr = fdr/100))
  assign(paste("candidates.",a,".",i,sep=""), which(get(paste("qv.",a,".",i,sep=""))$significant))
}}

# Saving candidates lists
for(a in c("lfmm","cate","emma","gapit")){
for(b in c("fisher","stouffer","weightedZ")){
  tmp=(get(paste("pv.",a,sep=""))[get(paste("candidates.",a,".",b,sep="")),c("SNP","Chr","Pos",paste("Combined_test_",b,sep=""))])
  names(tmp) <- c("SNP","Chr","Pos",paste0("pv.combined.",paste(pheno_variables,collapse = "_"),".",b))
  write.table(tmp,
                file = paste0("Combined_tests/candidates_",
                              paste(pheno_variables,collapse = "_"),"_",a,"_",b,"_fdr",fdr,"perc.txt"),
              quote = FALSE,
              row.names = FALSE)
}}

# Saving complete lists only for fisher results
for(a in c("lfmm","cate","emma","gapit")){
  for(b in c("fisher","stouffer","weightedZ")){
    tmp=(get(paste("pv.",a,sep=""))[,c("SNP","Chr","Pos",paste("Combined_test_",b,sep=""))])
    tmp <- cbind(tmp,get(paste0("qv.",a,".",b))$significant)
    names(tmp) <- c("SNP","Chr","Pos",paste0("pv.combined.",paste(pheno_variables,collapse = "_"),".",b),
                    paste0("signif.combined.",paste(pheno_variables,collapse = "_"),".",b,".",a,".fdr",fdr))
    write.table(tmp,
                file = paste0("Combined_tests/AllSNPs_",
                              paste(pheno_variables,collapse = "_"),"_",a,"_",b,"_fdr",fdr,"perc.txt"),
                quote = FALSE,
                row.names = FALSE)
  }
  
}

# Manhattan plots
for(a in c("lfmm","cate","emma","gapit")){
  CairoPNG(filename = paste("Combined_tests/man_plot_",a,"_",paste(pheno_variables,collapse = "_"),".png",sep=""),units = "cm",width = 21,height = 29.7,res = 600)
par(mfrow=c(3,2))
#for(i in pheno_variables){
#manhattan(x=get(paste("pv.",a,sep="")),chr = "Chromosome.gapit",bp = "Pos",p = paste("pv.",a,i,sep=""),
#          genomewideline = ifelse(test = length(get(paste("candidates.",a,".",i,sep="")))>0,
#                                  yes = -log10(max(get(paste("pv.",a,sep=""))[get(paste("candidates.",
#                                                                                        a,".",i,sep="")),paste("pv.",a,i,sep="")])),
#                                  no = FALSE),
#          suggestiveline = FALSE,
#          col=cols,main = paste(i))
#qq(pvector = unlist(get(paste("pv.",a,sep=""))[paste("pv.",a,i,sep="")]),main=paste("qqplot",a,i))
#  }
for(b in c("fisher","stouffer","weightedZ")){
manhattan(x=get(paste("pv.",a,sep="")),chr = "Chromosome.gapit",bp = "Pos",p = paste("Combined_test_",b,sep=""),
          genomewideline = ifelse(test = length(get(paste("candidates.",a,".",b,sep="")))>0,
                                  yes = -log10(max(get(paste("pv.",a,sep=""))[get(paste("candidates.",
                                                                                        a,".",b,sep="")),paste("Combined_test_",b,sep="")])),
                                  no = FALSE),
          suggestiveline = FALSE,
          col=cols,main = paste(b))
qq(pvector = unlist(get(paste("pv.",a,sep=""))[paste("Combined_test_",b,sep="")]),main=paste("qqplot",a,b,paste(pheno_variables,collapse = "_")))
}

dev.off()
}

for(a in c("lfmm","cate","emma","gapit")){
  pdf(file = paste("Combined_tests/man_plot_",a,"_",paste(pheno_variables,collapse = "_"),".pdf",sep=""),paper='a4')
  par(mfrow=c(3,2))
#  for(i in pheno_variables){
#    manhattan(x=get(paste("pv.",a,sep="")),chr = "Chromosome.gapit",bp = "Pos",p = paste("pv.",a,i,sep=""),
#              genomewideline = ifelse(test = length(get(paste("candidates.",a,".",i,sep="")))>0,
#                                      yes = -log10(max(get(paste("pv.",a,sep=""))[get(paste("candidates.",
#                                                                                            a,".",i,sep="")),paste("pv.",a,i,sep="")])),
#                                      no = FALSE),
#              suggestiveline = FALSE,
#              col=cols,main = paste(i))
#    qq(pvector = unlist(get(paste("pv.",a,sep=""))[paste("pv.",a,i,sep="")]),main=paste("qqplot",a,i))
#    
#  }
  for(b in c("fisher","stouffer","weightedZ")){
    manhattan(x=get(paste("pv.",a,sep="")),chr = "Chromosome.gapit",bp = "Pos",p = paste("Combined_test_",b,sep=""),
              genomewideline = ifelse(test = length(get(paste("candidates.",a,".",b,sep="")))>0,
                                      yes = -log10(max(get(paste("pv.",a,sep=""))[get(paste("candidates.",
                                                                                            a,".",b,sep="")),paste("Combined_test_",b,sep="")])),
                                      no = FALSE),
              suggestiveline = FALSE,
              col=cols,main = paste(b))
    qq(pvector = unlist(get(paste("pv.",a,sep=""))[paste("Combined_test_",b,sep="")]),main=paste("qqplot",a,b,paste(pheno_variables,collapse = "_")))
    
  }
  
  dev.off()
}
