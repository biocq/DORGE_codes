#####################
##### Bar plots #####
#####################

###### The boxplots to show the recall of CGC genes by prediction based on CRISPR-screening data only.

###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/TSG_prediction-lr_enet.CRISPR-only_FPR_0.01.csv: Predicted TSGs using CRISPR-screening data only under FPR = 0.01
###### Input: data/OG_prediction-lr_enet.CRISPR-only_FPR_0.01.csv: Predicted OGs using CRISPR-screening data only under FPR = 0.01
###### Output: Figure_2F_prediction_using_CRISPR_data_only.pdf: The number of recovered CGC-TSGs and CGC-OGs using all features compared to CRISPR-screening data only

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggplot2")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("ggplot2"))


TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

################################### Figure 2F ###################################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_2");
CRISPR_specific <- read.table("data/TSG_prediction-lr_enet.CRISPR-only_FPR_0.01.csv", header=T, sep=",")
Predicted_genes_CRISPR_features<-as.character(CRISPR_specific[,1])
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
index_TSG_predicted_CGC_TSGs_all_feature<-which(prediction$TSG_probability>TSG_threshold & prediction$TSG_all=="1")
all_genes<-as.character(prediction[,1])
Predicted_CGC_TSGs_all_features<-all_genes[index_TSG_predicted_CGC_TSGs_all_feature]
index_predicted_CRISPR_genes_in_all_genes<-which(all_genes %in% Predicted_genes_CRISPR_features)
index_CGC_TSG<-which(prediction$TSG_all==1)
Predicted_CGC_TSGs_CRISPR_features_only<-all_genes[intersect(index_predicted_CRISPR_genes_in_all_genes,index_CGC_TSG)]
#CRISPR_specific_genes<-setdiff(Predicted_CGC_TSGs_CRISPR_features_only,Predicted_CGC_TSGs_all_features)

CRISPR_specific <- read.table("data/OG_prediction-lr_enet.CRISPR-only_FPR_0.01.csv", header=T, sep=",")
Predicted_genes_CRISPR_features<-as.character(CRISPR_specific[,1])
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
index_OG_predicted_CGC_OGs_all_feature<-which(prediction$OG_probability>OG_threshold & prediction$OG_all=="1")
all_genes<-as.character(prediction[,1])
Predicted_CGC_OGs_all_features<-all_genes[index_OG_predicted_CGC_OGs_all_feature]
index_predicted_CRISPR_genes_in_all_genes<-which(all_genes %in% Predicted_genes_CRISPR_features)
index_CGC_OG<-which(prediction$OG_all==1)
Predicted_CGC_OGs_CRISPR_features_only<-all_genes[intersect(index_predicted_CRISPR_genes_in_all_genes,index_CGC_OG)]
#CRISPR_specific_genes<-setdiff(Predicted_CGC_OGs_CRISPR_features_only,Predicted_CGC_OGs_all_features)

dat <- data.frame("Feature" = c(length(Predicted_CGC_TSGs_all_features)-length(Predicted_CGC_TSGs_CRISPR_features_only),length(Predicted_CGC_TSGs_CRISPR_features_only),length(Predicted_CGC_OGs_all_features)-length(Predicted_CGC_OGs_CRISPR_features_only),length(Predicted_CGC_OGs_CRISPR_features_only)), "Group" = c("All-feature-specific","CRISPR-screening-only","All-feature-specific","CRISPR-screening-only"),"Cat" = c("CGC-TSG","CGC-TSG","CGC-OG","CGC-OG"))
dat$Cat<-factor(dat$Cat,levels=c("CGC-TSG", "CGC-OG"))

pdf("Raw_figures/Figure_2F_prediction_using_CRISPR_data_only.pdf", family="ArialMT", width=3.5, height=3)
ggplot(dat, aes(fill=Group, y=Feature, x=Cat)) + geom_bar(position="stack", stat="identity")+ labs(x = "", y = "Gene number") +scale_fill_manual(values=c("red", "purple"))+ theme_bw() + theme(legend.title = element_blank(),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
garbage <- dev.off()