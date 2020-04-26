###########################################
###### Prediction result processing  ######
###########################################

###### The feature profile "All_features.csv" is used as the annotation sources.
###### Input: ../../Gene_set_new.txt: Gene annotation file compiled in this study
###### Input: TSG_prediction-lr_enet.csv: TSG prediction results by elastic net
###### Input: OG_prediction-lr_enet.csv: OG prediction results by elastic net
###### Output: ../../DORGE_prediction.txt: Formatted DORGE prediction table including both TSG probability and OG probability as well as annotation from Gene_set_new.txt
###### Output: ../Novel_DORGE_predicted_TSGs.txt: Gene list of novel non-CGC DORGE-predicted TSGs.
###### Output: ../Novel_DORGE_predicted_OGs.txt: Gene list of novel non-CGC DORGE-predicted OGs.

### Install missing packages
installed_pkgs <- installed.packages()
pkgs <-  c("ComplexHeatmap")

pkgs <-  c("plyr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

suppressMessages(library(plyr))

options(warn=-1)
################################### Data file S2: Formatted prediction results ###################################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/Data_file_S2/");

anno <- read.table("../../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
TSG_result <- read.table("TSG_prediction-lr_enet.csv", header=T, sep=",",fill=TRUE,quote = "")
OG_result <- read.table("OG_prediction-lr_enet.csv", header=T, sep=",",fill=TRUE,quote = "")
join<-join(anno, TSG_result, type = "left", by="Gene")
join<-join(join, OG_result, type = "left", by="Gene")
write.table(join[,c(1,7,8,2,3,4,5,6)],"../../DORGE_prediction.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(join[,c(1,7,8,2,3,4,5,6)],"DORGE_prediction.txt",sep="\t",quote=F,row.names=F,col.names=T)

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]

index_TSG<-which(TSG_CGC=="1")
index_OG<-which(OG_CGC=="1")
index_pTSG<-which((join$TSG_probability>0.62485 & TSG_CGC!="1"))#FPR<0.01
index_pOG<-which((join$OG_probability>0.7004394 & OG_CGC!="1"))#FPR<0.01


Novel_TSG<-as.character(join[index_pTSG,1])
Novel_OG<-as.character(join[index_pOG,1])
write.table(Novel_TSG,"../Novel_DORGE_predicted_TSGs.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Novel_OG,"../Novel_DORGE_predicted_OGs.txt",sep="\t",quote=F,row.names=F,col.names=F)

#write.table(join[index_pTSG,1:8],"../Novel_DORGE_predicted_TSGs1.txt",sep="\t",quote=F,row.names=F,col.names=T)
#write.table(join[index_pOG,1:8],"../Novel_DORGE_predicted_OGs1.txt",sep="\t",quote=F,row.names=F,col.names=T)