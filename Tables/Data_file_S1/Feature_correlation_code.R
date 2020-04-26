##############################################
###### Feature correlation calculation  ######
##############################################

###### Input: ../../All_features.csv: Feature table compiled in this study
###### Output: corr_matrix.txt: Feature correlation table

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("Hmisc")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs),repos = "http://cran.us.r-project.org")
}

suppressMessages(library(Hmisc))
################# Correlation between pair-wise features #################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/Data_file_S1/");

all_feature <- read.table("../../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
subset<-c(2:23,25:33,35:37,39:40,42,44,48,50,52:53,55:66,68,69,71,72,74,75,77,78,80,81,83,84,86,87,89,90,92,93,95,96,98)
res2 <- rcorr(as.matrix(all_feature[,subset]))
m3<-res2$r
diag(m3) <- 0
write.table(m3,"corr_matrix.txt",sep="\t",quote=F,row.names=T)
#m4<-res2$P
#write.table(m4,"corr_matrix_p.txt",sep="\t",quote=F,row.names=T)