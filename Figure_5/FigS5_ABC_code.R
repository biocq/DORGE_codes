####################
##### Boxplots #####
####################

###### The boxplots for different gene categories to independently evaluate DORGE prediction.

###### Input: ../Gene_set_new.txt: Gene annotation file
###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/BioGRID_network_node_metrics_processed.txt:  Network metrics calculated by Cytoscape (https://cytoscape.org/; Available at Evaluation_processing folder DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)

###### Output: Figure_S5A_BioGRID_log_degree.pdf: Evaluation of dual-function TSG/OGs by PPI network degree
###### Output: Figure_S5B_BioGRID_betweenness.pdf: Evaluation of dual-function TSG/OGs by PPI network betweenness centrality
###### Output: Figure_S5C_BioGRID_closeness.pdf: Evaluation of dual-function TSG/OGs by PPI network closeness centrality

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggsignif","ggpubr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

suppressMessages(library("ggpubr"))
suppressMessages(library(ggsignif))

TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

######### Figure S6 ABC: Network metrics in PPI network #########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_5");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
anno<-cbind(anno[,c(1,4,5,6)],prediction[,2:3])
all_feature <- read.table("data/BioGRID_network_node_metrics_processed.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","Gene_betweenness_centrality","Gene_closeness_centrality","Gene_Degree")]
all_feature$Gene_betweenness_centrality <- as.numeric(all_feature$Gene_betweenness_centrality)
all_feature$Gene_betweenness_centrality <- log10(all_feature$Gene_betweenness_centrality+1e-8)
all_feature$Gene_closeness_centrality <- as.numeric(all_feature$Gene_closeness_centrality)
all_feature$Gene_Degree <- as.numeric(all_feature$Gene_Degree)
all_feature$Gene_Degree <- log2(all_feature$Gene_Degree+1)

colnames(all_feature)<-c("Gene","Log10_BioGRID_betweenness","BioGRID_clossness","Log_BioGRID_degree")
index_low<-which(all_feature$Log10_BioGRID_betweenness< -6)
all_feature$Log10_BioGRID_betweenness[index_low]<-(all_feature$Log10_BioGRID_betweenness[index_low]+6)/10-6

index_low<-which(all_feature$BioGRID_clossness<0.3)
all_feature$BioGRID_clossness[index_low]<-(all_feature$BioGRID_clossness[index_low]-0.3)/10+0.3
index_high<-which(all_feature$BioGRID_clossness>0.5)
all_feature$BioGRID_clossness[index_high]<-(all_feature$BioGRID_clossness[index_high]-0.5)/100+0.5

index_high<-which(all_feature$Log_BioGRID_degree>7)
all_feature$Log_BioGRID_degree[index_high]<-(all_feature$Log_BioGRID_degree[index_high]-7)/5+7

index<-rep("Other",nrow(anno))
index_NG<-which(anno$NG=="1")
index_pOG<-which(anno$OG_probability> OG_threshold & anno$TSG_probability<TSG_threshold & anno$OG_all!="1")
index_pTSG<-which(anno$TSG_probability>TSG_threshold & anno$OG_probability< OG_threshold & anno$TSG_all!="1")
index_pdual<-which((anno$OG_all=="1" & anno$TSG_all!="1" & anno$TSG_probability>TSG_threshold)|(anno$TSG_all=="1" & anno$OG_all!="1" & anno$OG_probability> OG_threshold)|(anno$TSG_all!="1" & anno$OG_all!="1" & anno$TSG_probability>TSG_threshold & anno$OG_probability> OG_threshold))
index_TSG<-which(anno$TSG_all=="1" & anno$OG_all!="1")
index_OG<-which(anno$OG_all=="1" & anno$TSG_all!="1")
index_dual<-which(anno$OG_all=="1" & anno$TSG_all=="1")

index[index_NG]<-"NG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"
index[index_pdual]<-"Novel DORGE-dual"
index[index_dual]<-"CGC-dual"

dat2<-cbind(all_feature,index)
dat3<-dat2[dat2$index!="Other",]
dat3$index<-factor(dat3$index, levels=c("CGC-dual","NG","CGC-OG","CGC-TSG","Novel DORGE-OG","Novel DORGE-TSG","Novel DORGE-dual"))

pdf("Raw_figures/Figure_S5A_BioGRID_log_degree.pdf", family="ArialMT", width=2.85, height=5) #Degree
p<-ggboxplot(data = dat3,x = "index", y="Log_BioGRID_degree",color="black",fill="index", palette = rev(c("#FF00FF","#984EA3","#4DAF4A","#377EB8","#E41A1C","#CCCCCC","#999999")),outlier.size=0.2,width=0.7,lwd=0.25)+scale_y_continuous(breaks = c(0,3,7),labels = c("0","3","7"))+ labs(x = "", y = "log Degree") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.y = element_text(hjust = 1, colour = "black"),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="NG","Log_BioGRID_degree"]))+0.8, max(max(dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="CGC-OG","Log_BioGRID_degree"]))+1.3),xmin=c(1.05,1.05), xmax=c(1.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"],dat3[dat3$index=="NG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"],dat3[dat3$index=="CGC-OG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"]))+0.3, max(max(dat3[dat3$index=="CGC-OG","Log_BioGRID_degree"]),max(dat3[dat3$index=="NG","Log_BioGRID_degree"]))+0.1),xmin=c(1.05,2.05), xmax=c(3.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"],dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-OG","Log_BioGRID_degree"],dat3[dat3$index=="NG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="Novel DORGE-OG","Log_BioGRID_degree"]))+1.2, max(max(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="Novel DORGE-TSG","Log_BioGRID_degree"]))+0.3),xmin=c(5.05,6.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"],dat3[dat3$index=="Novel DORGE-OG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"],dat3[dat3$index=="Novel DORGE-TSG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","Log_BioGRID_degree"]),max(dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"]))+0.6, max(max(dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"]),max(dat3[dat3$index=="NG","Log_BioGRID_degree"]))+2.0),xmin=c(4.05,2.05), xmax=c(4.95,3.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","Log_BioGRID_degree"],dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"],dat3[dat3$index=="NG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)

p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-TSG","Log_BioGRID_degree"]),max(dat3[dat3$index=="CGC-OG","Log_BioGRID_degree"]))+2.7, max(max(dat3[dat3$index=="Novel DORGE-TSG","Log_BioGRID_degree"]),max(dat3[dat3$index=="NG","Log_BioGRID_degree"]))+3.2),xmin=c(3.05,2.05), xmax=c(5.95,5.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","Log_BioGRID_degree"],dat3[dat3$index=="CGC-OG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","Log_BioGRID_degree"],dat3[dat3$index=="NG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","Log_BioGRID_degree"]),max(dat3[dat3$index=="NG","Log_BioGRID_degree"]))+4.2, max(max(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="NG","Log_BioGRID_degree"]))+3.6),xmin=c(2.05,2.05), xmax=c(4.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","Log_BioGRID_degree"],dat3[dat3$index=="NG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"],dat3[dat3$index=="NG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"]))+4.6, max(max(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"]),max(dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"]))+1.8),xmin=c(1.05,4.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"],dat3[dat3$index=="CGC-dual","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log_BioGRID_degree"],dat3[dat3$index=="CGC-TSG","Log_BioGRID_degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p
garbage <- dev.off()

# Column: 1 CGCdual  2 NG  3 CGCOG  4 CGCTSG  5 Novel DORGEOG  6 Novel DORGETSG  7 Novel DORGEdual
#CGC-dual NG; CGC-dual CGC-OG
#CGC-dual CGC-TSG; CGC-OG NG
#Novel DORGE-dual Novel DORGE-OG; Novel DORGE-dual Novel DORGE-TSG
#Novel DORGE-OG CGC-TSG; CGC-TSG NG
#Novel DORGE-TSG CGC-OG; Novel DORGE-TSG NG
#Novel DORGE-OG NG; Novel DORGE-dual NG
#Novel DORGE-dual CGC-dual; Novel DORGE-dual CGC-TSG

pdf("Raw_figures/Figure_S5B_BioGRID_betweenness.pdf", family="ArialMT", width=2.85, height=5) #Betweenness Centrality
p<-ggboxplot(data = dat3,x = "index", y="Log10_BioGRID_betweenness",color="black",fill="index", palette = rev(c("#FF00FF","#984EA3","#4DAF4A","#377EB8","#E41A1C","#CCCCCC","#999999")),outlier.size=0.2,width=0.7,lwd=0.25)+scale_y_continuous(breaks = c(-6,-3,-1),labels = c("-6","-3","-1"))+ labs(x = "", y = "log10(Betweenness centrality)") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.y = element_text(hjust = 1, colour = "black"),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="NG","Log10_BioGRID_betweenness"]))+0.8/2, max(max(dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="CGC-OG","Log10_BioGRID_betweenness"]))+1.3/2),xmin=c(1.05,1.05), xmax=c(1.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="NG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="CGC-OG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"]))+0.3/2, max(max(dat3[dat3$index=="CGC-OG","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="NG","Log10_BioGRID_betweenness"]))+0.1/2),xmin=c(1.05,2.05), xmax=c(3.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-OG","Log10_BioGRID_betweenness"],dat3[dat3$index=="NG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="Novel DORGE-OG","Log10_BioGRID_betweenness"]))+0.8/2, max(max(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="Novel DORGE-TSG","Log10_BioGRID_betweenness"]))+0.3/2),xmin=c(5.05,6.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="Novel DORGE-OG","BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="Novel DORGE-TSG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"]))+0.6/2, max(max(dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="NG","Log10_BioGRID_betweenness"]))+2.6/2),xmin=c(4.05,2.05), xmax=c(4.95,3.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","Log10_BioGRID_betweenness"],dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"],dat3[dat3$index=="NG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)

p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-TSG","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="CGC-OG","Log10_BioGRID_betweenness"]))+2.7/2, max(max(dat3[dat3$index=="Novel DORGE-TSG","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="NG","Log10_BioGRID_betweenness"]))+3.2/2),xmin=c(3.05,2.05), xmax=c(5.95,5.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","Log10_BioGRID_betweenness"],dat3[dat3$index=="CGC-OG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","Log10_BioGRID_betweenness"],dat3[dat3$index=="NG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="NG","Log10_BioGRID_betweenness"]))+3.6/2, max(max(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="NG","Log10_BioGRID_betweenness"]))+3.4/2),xmin=c(2.05,2.05), xmax=c(4.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","Log10_BioGRID_betweenness"],dat3[dat3$index=="NG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="NG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"]))+4.6/2, max(max(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"]),max(dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"]))+1.8/2),xmin=c(1.05,4.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="CGC-dual","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Log10_BioGRID_betweenness"],dat3[dat3$index=="CGC-TSG","Log10_BioGRID_betweenness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p
garbage <- dev.off()

pdf("Raw_figures/Figure_S5C_BioGRID_closeness.pdf", family="ArialMT", width=2.85, height=5) #Closeness Centrality
p<-ggboxplot(data = dat3,x = "index", y="BioGRID_clossness",color="black",fill="index", palette = rev(c("#FF00FF","#984EA3","#4DAF4A","#377EB8","#E41A1C","#CCCCCC","#999999")),outlier.size=0.2,width=0.7,lwd=0.25)+scale_y_continuous(breaks = c(0.3,0.4,0.5),labels = c("0.3","0.4","0.5"))+ labs(x = "", y = "Closeness centrality") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.y = element_text(hjust = 1, colour = "black"),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","BioGRID_clossness"]),max(dat3[dat3$index=="NG","BioGRID_clossness"]))+0.8/25, max(max(dat3[dat3$index=="CGC-dual","BioGRID_clossness"]),max(dat3[dat3$index=="CGC-OG","BioGRID_clossness"]))+1.3/25),xmin=c(1.05,1.05), xmax=c(1.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","BioGRID_clossness"],dat3[dat3$index=="NG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-dual","BioGRID_clossness"],dat3[dat3$index=="CGC-OG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","BioGRID_clossness"]),max(dat3[dat3$index=="CGC-TSG","BioGRID_clossness"]))+0.3/25, max(max(dat3[dat3$index=="CGC-OG","BioGRID_clossness"]),max(dat3[dat3$index=="NG","BioGRID_clossness"]))+0.4/25),xmin=c(1.05,2.05), xmax=c(3.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","BioGRID_clossness"],dat3[dat3$index=="CGC-TSG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-OG","BioGRID_clossness"],dat3[dat3$index=="NG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"]),max(dat3[dat3$index=="Novel DORGE-OG","BioGRID_clossness"]))+1.2/25, max(max(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"]),max(dat3[dat3$index=="Novel DORGE-TSG","BioGRID_clossness"]))+0.3/25),xmin=c(5.05,6.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"],dat3[dat3$index=="Novel DORGE-OG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"],dat3[dat3$index=="Novel DORGE-TSG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","BioGRID_clossness"]),max(dat3[dat3$index=="CGC-TSG","BioGRID_clossness"]))+0.6/25, max(max(dat3[dat3$index=="CGC-TSG","BioGRID_clossness"]),max(dat3[dat3$index=="NG","BioGRID_clossness"]))+2.0/25),xmin=c(4.05,2.05), xmax=c(4.95,3.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","BioGRID_clossness"],dat3[dat3$index=="CGC-TSG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-TSG","BioGRID_clossness"],dat3[dat3$index=="NG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)

p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-TSG","BioGRID_clossness"]),max(dat3[dat3$index=="CGC-OG","BioGRID_clossness"]))+2.7/25, max(max(dat3[dat3$index=="Novel DORGE-TSG","BioGRID_clossness"]),max(dat3[dat3$index=="NG","BioGRID_clossness"]))+3.2/25),xmin=c(3.05,2.05), xmax=c(5.95,5.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","BioGRID_clossness"],dat3[dat3$index=="CGC-OG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","BioGRID_clossness"],dat3[dat3$index=="NG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","BioGRID_clossness"]),max(dat3[dat3$index=="NG","BioGRID_clossness"]))+4.2/25, max(max(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"]),max(dat3[dat3$index=="NG","BioGRID_clossness"]))+3.6/25),xmin=c(2.05,2.05), xmax=c(4.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","BioGRID_clossness"],dat3[dat3$index=="NG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"],dat3[dat3$index=="NG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"]),max(dat3[dat3$index=="CGC-dual","BioGRID_clossness"]))+5.0/25, max(max(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"]),max(dat3[dat3$index=="CGC-TSG","BioGRID_clossness"]))+1.8/25),xmin=c(1.05,4.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"],dat3[dat3$index=="CGC-dual","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","BioGRID_clossness"],dat3[dat3$index=="CGC-TSG","BioGRID_clossness"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p
garbage <- dev.off()