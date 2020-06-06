###########################################
###### Prediction result processing  ######
###########################################

###### Input: ../../Gene_set_new.txt: Gene annotation file compiled in this study
###### Input: TSG_prediction-lr_enet.wo_Epigenetics.csv: TSG prediction results predicted by DORGE based on non-Epigenetic features.
###### Input: OG_prediction-lr_enet.wo_Epigenetics.csv: OG prediction results predicted by DORGE based on non-Epigenetic features.
###### Output: Novel_DORGE_predicted_TSGs_wo_Epi.txt: Gene list of novel non-CGC DORGE-predicted TSGs predicted by DORGE based on non-Epigenetic features.
###### Output: Novel_DORGE_predicted_OGs_wo_Epi.txt: Gene list of novel non-CGC DORGE-predicted OGs predicted by DORGE based on non-Epigenetic features.

### Install missing packages
installed_pkgs <- installed.packages()

pkgs <-  c("plyr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library(plyr))

options(warn=-1)

TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

################################### Data file S2: Formatted prediction results ###################################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_3/data/");

anno <- read.table("../../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
TSG_result <- read.table("TSG_prediction-lr_enet.wo_Epigenetics.csv", header=T, sep=",",fill=TRUE,quote = "")
OG_result <- read.table("OG_prediction-lr_enet.wo_Epigenetics.csv", header=T, sep=",",fill=TRUE,quote = "")
join<-join(anno, TSG_result, type = "left", by="Gene")
join<-join(join, OG_result, type = "left", by="Gene")
#write.table(join[,c(1,7,8,2,3,4,5,6)],"../../DORGE_prediction.txt",sep="\t",quote=F,row.names=F,col.names=T)

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]

index_TSG<-which(TSG_CGC=="1")
index_OG<-which(OG_CGC=="1")
index_pTSG<-which((join$TSG_probability>TSG_threshold & TSG_CGC!="1"))#FPR<0.01
index_pOG<-which((join$OG_probability>OG_threshold & OG_CGC!="1"))#FPR<0.01


Novel_TSG<-as.character(join[index_pTSG,1])
Novel_OG<-as.character(join[index_pOG,1])
write.table(Novel_TSG,"Novel_DORGE_predicted_TSGs_wo_Epi.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Novel_OG,"Novel_DORGE_predicted_OGs_wo_Epi.txt",sep="\t",quote=F,row.names=F,col.names=F)

