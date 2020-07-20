####################
##### Boxplots #####
####################

###### The boxplots for survival data to independently evaluate DORGE prediction.

###### Input: ../Gene_set_new.txt: Gene annotation file
###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/oncolnc_processed.txt: Processed survival data (from https://peerj.com/articles/cs-67/; Available at Evaluation_processing folder DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)

###### Output: Figure_3E_survival_READ_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by survival data in READ
###### Output: Figure_S4E_survival_COAD_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by survival data in COAD
###### Output: Figure_S4F_survival_UCEC_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by survival data in UCEC

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("plyr","ggpubr")


if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("ggpubr"))
suppressMessages(library("plyr"))

TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

##################### Figure 3E, S4E and S4F: Survival data evaluation #####################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_3");

anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
anno<-cbind(anno[,c(1,4,5,6)],prediction[,2:3])
all_feature <- read.table("data/oncolnc_all_genes.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","cancer","HR")]
all_feature$HR <- as.numeric(all_feature$HR)

join<-join(all_feature, anno,type = "full",by="Gene")
join<-join[!is.na(join$HR),]
index_pOG<-which(join$OG_probability> OG_threshold & join$TSG_probability< TSG_threshold & join$OG_all!="1")
index_pTSG<-which(join$TSG_probability> TSG_threshold & join$OG_probability< OG_threshold & join$TSG_all!="1")
index_TSG<-which(join$TSG_all=="1" & join$OG_all!="1")
index_OG<-which(join$OG_all=="1" & join$TSG_all!="1")
index_NG<-which(join$NG=="1")

index<-rep("",nrow(join))

index[index_NG]<-"NG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

join<-cbind(join,index)
allgene2<-join[join$index!="" & !is.na(join$index),]
allgene2$index<-factor(allgene2$index,levels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))

ca<-paste("Raw_figures/Figure_3E_survival_READ_evaluation.pdf",sep="")
pdf(ca, family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-OG","CGC-TSG"),c("Novel DORGE-OG","Novel DORGE-TSG"),c("NG","CGC-TSG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("NG","Novel DORGE-TSG"),c("Novel DORGE-OG","CGC-TSG"),c("CGC-OG","Novel DORGE-TSG"));
ggboxplot(data = allgene2[allgene2$cancer=="READ",],x = "index", y="HR",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_y_continuous(breaks = c(0.5,1.0,1.5),labels = c("0.5", "1.0", "1.5"))+labs(x = "", y = "Hazard ratio")+ggtitle("READ") + theme_bw() + theme(plot.title = element_text(size = 5),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

ca<-paste("Raw_figures/Figure_S4E_survival_COAD_evaluation.pdf",sep="")
pdf(ca, family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-OG","CGC-TSG"),c("Novel DORGE-OG","Novel DORGE-TSG"),c("NG","CGC-TSG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("NG","Novel DORGE-TSG"),c("Novel DORGE-OG","CGC-TSG"),c("CGC-OG","Novel DORGE-TSG"));
ggboxplot(data = allgene2[allgene2$cancer=="COAD",],x = "index", y="HR",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_y_continuous(breaks = c(0.5,1.0,1.5),labels = c("0.5", "1.0", "1.5"))+labs(x = "", y = "Hazard ratio")+ggtitle("COAD") + theme_bw() + theme(plot.title = element_text(size = 5),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

ca<-paste("Raw_figures/Figure_S4F_survival_UCEC_evaluation.pdf",sep="")
pdf(ca, family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-OG","CGC-TSG"),c("Novel DORGE-OG","Novel DORGE-TSG"),c("NG","CGC-TSG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("NG","Novel DORGE-TSG"),c("Novel DORGE-OG","CGC-TSG"),c("CGC-OG","Novel DORGE-TSG"));
ggboxplot(data = allgene2[allgene2$cancer=="UCEC",],x = "index", y="HR",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_y_continuous(breaks = c(0.5,1.0,1.5),labels = c("0.5", "1.0", "1.5"))+labs(x = "", y = "Hazard ratio")+ggtitle("UCEC") + theme_bw() + theme(plot.title = element_text(size = 5),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()



cancer_types<-as.character(unique(allgene2$cancer))
for(i in cancer_types){
	#ca<-paste("Raw_figures/Figure_S4E_oncolnc_",i,"_evaluation.pdf",sep="")
	#pdf(ca, family="ArialMT", width=2.275, height=3)
	#my_comparisons <- list(c("CGC-OG","CGC-TSG"),c("Novel DORGE-OG","Novel DORGE-TSG"),c("NG","CGC-TSG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("NG","Novel DORGE-TSG"),c("Novel DORGE-OG","CGC-TSG"),c("CGC-OG","Novel DORGE-TSG"));
	#print(ggboxplot(data = allgene2[allgene2$cancer==i,],x = "index", y="HR",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_y_continuous(breaks = c(0.5,1.0,1.5),labels = c("0.5", "1.0", "1.5"))+labs(x = "", y = "Hazard ratio")+ggtitle(i) + theme_bw() + theme(plot.title = element_text(size = 5),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5))
	#garbage <- dev.off()
}