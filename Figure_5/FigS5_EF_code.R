###############################
##### Boxplot and barplot #####
###############################

###### The boxplots for different gene categories on PharmacoDB network to independently evaluate DORGE prediction.

###### Input: ../Gene_set_new.txt: Gene annotation file
###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/gene_degree_pharmacodb.txt: Degree metric calculated based on PharmacoDB data (Available at Evaluation_processing folder DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)

###### Output: Figure_S5E_pharmacodb_log_degree.pdf: Evaluation of CGC and DORGE-predicted dual-function TSG/OGs by Pharmacodb network degree
###### Output: Figure_S5F_drug_degree.pdf: Top 10 drugs/compounds that are most associated with genes in Pharmacodb

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("plyr","ggsignif","ggpubr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library(ggpubr))
suppressMessages(library(ggsignif))
suppressMessages(library(plyr))

TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

######### Figure S5E: Network metrics in drug-compound network #########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_5");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
anno<-cbind(anno[,c(1,4,5,6)],prediction[,2:3])
all_feature <- read.table("data/gene_degree_pharmacodb.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","Degree")]
all_feature$Degree <- as.numeric(all_feature$Degree)
all_feature$Degree <- log2(all_feature$Degree+1)
index_high<-which(all_feature$Degree>6)
all_feature$Degree[index_high]<-(all_feature$Degree[index_high]-6)/10+6

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

pdf("Raw_figures/Figure_S5E_pharmacodb_log_degree.pdf", family="ArialMT", width=2.85, height=5) #Degree
p<-ggboxplot(data = dat3,x = "index", y="Degree",color="black",fill="index", palette = rev(c("#FF00FF","#984EA3","#4DAF4A","#377EB8","#E41A1C","#CCCCCC","#999999")),outlier.size=0.2,width=0.7,lwd=0.25)+scale_y_continuous(breaks = c(0,3,7),labels = c("0","3","7"))+ labs(x = "", y = "log Degree") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.y = element_text(hjust = 1, colour = "black"),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","Degree"]),max(dat3[dat3$index=="NG","Degree"]))+0.8, max(max(dat3[dat3$index=="CGC-dual","Degree"]),max(dat3[dat3$index=="CGC-OG","Degree"]))+1.3),xmin=c(1.05,1.05), xmax=c(1.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Degree"],dat3[dat3$index=="NG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Degree"],dat3[dat3$index=="CGC-OG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="CGC-dual","Degree"]),max(dat3[dat3$index=="CGC-TSG","Degree"]))+0.5, max(max(dat3[dat3$index=="CGC-OG","Degree"]),max(dat3[dat3$index=="NG","Degree"]))+0.1),xmin=c(1.05,2.05), xmax=c(3.95,2.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="CGC-dual","Degree"],dat3[dat3$index=="CGC-TSG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-OG","Degree"],dat3[dat3$index=="NG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","Degree"]),max(dat3[dat3$index=="Novel DORGE-OG","Degree"]))+1.2, max(max(dat3[dat3$index=="Novel DORGE-dual","Degree"]),max(dat3[dat3$index=="Novel DORGE-TSG","Degree"]))+0.3),xmin=c(5.05,6.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Degree"],dat3[dat3$index=="Novel DORGE-OG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Degree"],dat3[dat3$index=="Novel DORGE-TSG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","Degree"]),max(dat3[dat3$index=="CGC-TSG","Degree"]))+0.6, max(max(dat3[dat3$index=="CGC-TSG","Degree"]),max(dat3[dat3$index=="NG","Degree"]))+2.0),xmin=c(4.05,2.05), xmax=c(4.95,3.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","Degree"],dat3[dat3$index=="CGC-TSG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="CGC-TSG","Degree"],dat3[dat3$index=="NG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)

p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-TSG","Degree"]),max(dat3[dat3$index=="CGC-OG","Degree"]))+2.7, max(max(dat3[dat3$index=="Novel DORGE-TSG","Degree"]),max(dat3[dat3$index=="NG","Degree"]))+3.2),xmin=c(3.05,2.05), xmax=c(5.95,5.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","Degree"],dat3[dat3$index=="CGC-OG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-TSG","Degree"],dat3[dat3$index=="NG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-OG","Degree"]),max(dat3[dat3$index=="NG","Degree"]))+4.2, max(max(dat3[dat3$index=="Novel DORGE-dual","Degree"]),max(dat3[dat3$index=="NG","Degree"]))+3.6),xmin=c(2.05,2.05), xmax=c(4.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-OG","Degree"],dat3[dat3$index=="NG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Degree"],dat3[dat3$index=="NG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p<-p+geom_signif(y_position=c(max(max(dat3[dat3$index=="Novel DORGE-dual","Degree"]),max(dat3[dat3$index=="CGC-dual","Degree"]))+4.6, max(max(dat3[dat3$index=="Novel DORGE-dual","Degree"]),max(dat3[dat3$index=="CGC-TSG","Degree"]))+1.8),xmin=c(1.05,4.05), xmax=c(6.95,6.95),annotation=c(formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Degree"],dat3[dat3$index=="CGC-dual","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2), formatC(wilcox.test(dat3[dat3$index=="Novel DORGE-dual","Degree"],dat3[dat3$index=="CGC-TSG","Degree"],alternative = "two.sided")$p.value, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =3)
p
garbage <- dev.off()

######## Figure S5F Drug/compound degree #########


#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_5");
drug_gene_nr <- read.table("data/drug_gene_nr.txt", header=T, sep="\t")

count<-plyr::count(drug_gene_nr, "Drug")
count<-count[rev(order(as.numeric(count$freq))),]

pdf("Raw_figures/Figure_S5F_drug_degree.pdf", family="ArialMT", width=3, height=4)
ggbarplot(count[1:10,],x="Drug",y="freq",sort.val = "desc",fill = "steelblue",color="transparent",rotate = TRUE,width = 0.8)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"))+labs(x="",y="Gene number")+ guides(fill=FALSE)
garbage <- dev.off()