####################
##### Boxplots #####
####################

###### The boxplots for different gene categories to independently evaluate DORGE prediction.

###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: sampled_NGs.txt: Sampled 2,000 NGs for reproducing purpose
###### Input: data/pan_cancer_ATACseq_peaks_processed.txt: ATAC-seq pan-cancer data (Data S2 of https://science.sciencemag.org/content/362/6413/eaav1898; Available at Evaluation_processing folder from the DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)
###### Input: data/Noncoding_mutation_rate_processed.txt: noncoding mutation rate data (Processed based on Table S5 of https://www.nature.com/articles/nature12213; Available at Evaluation_processing folder from the DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)
###### Input: data/shRNA_cell_proliferation_rates.txt: DepMap shRNA screening data data (Achilles v2.4.5; Available at Evaluation_processing folder from the DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)
###### Input: data/phyloP46_processed.txt: phyloP data downloaded from UCSC (Available at Evaluation_processing folder from the DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)

###### Output: Figure_4A_ATAC_seq_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by ATAC-seq data
###### Output: Figure_S4F_noncoding_mutation_evaluation.pdf:  Evaluation of DORGE novel TSG/OGs by noncoding mutation rate data
###### Output: Figure_S4C_shRNA_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by shRNA screening data
###### Output: Figure_S4D_phastcons_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by phastcons feature
###### Output: Figure_4E_PhyloP_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by PhyloP feature
###### Output: Figure_S4E_ncRVIS_evaluation.pdf: Evaluation of DORGE novel TSG/OGs by ncRVIS feature

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggpubr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("ggpubr"))

TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

######### Figure 4A: ATAC-seq data evaluation for DORGE prediction #########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_4");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("data/pan_cancer_ATACseq_peaks_processed.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","ATAC_seq_peak_signal")]
all_feature$ATAC_seq_peak_signal <- as.numeric(all_feature$ATAC_seq_peak_signal)
index_large_ATAC_signal<-which(all_feature$ATAC_seq_peak_signal>200)
all_feature$ATAC_seq_peak_signal[index_large_ATAC_signal]<-(all_feature$ATAC_seq_peak_signal[index_large_ATAC_signal]-200)/25+200

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")

index_pOG<-which(prediction$OG_probability> OG_threshold & prediction$TSG_probability<TSG_threshold & OG_CGC!="1")
index_pTSG<-which(prediction$TSG_probability>TSG_threshold & prediction$OG_probability< OG_threshold & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")

index<-rep("",nrow(anno))
index[which(NG=="1")]<-"NG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

all_feature<-cbind(all_feature,index)
allgene2<-all_feature[all_feature$index!="",]
allgene2$index<-factor(allgene2$index,levels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))
pdf("Raw_figures/Figure_4A_ATAC_seq_evaluation.pdf", family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-OG","CGC-TSG"),c("Novel DORGE-OG","Novel DORGE-TSG"),c("CGC-TSG","NG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("Novel DORGE-TSG","NG"),c("Novel DORGE-OG","CGC-TSG"),c("CGC-OG","Novel DORGE-TSG"));

ggboxplot(data = allgene2,x = "index", y="ATAC_seq_peak_signal",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_x_discrete(labels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))+scale_y_continuous(breaks = c(0,100,200,300),labels = c("0", "100", "200", "2700"))+labs(x = "", y = "ATAC-seq peak score") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <-dev.off()

######## Figure S4F noncoding_mutation_rate evaluation########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_4");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("data/Noncoding_mutation_rate_processed.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","noncoding_mutation_rate")]
all_feature$noncoding_mutation_rate <- as.numeric(all_feature$noncoding_mutation_rate)
index_largenoncoding_mutation_rate<-which(all_feature$noncoding_mutation_rate> 4e-6)
index_lownoncoding_mutation_rate<-which(all_feature$noncoding_mutation_rate< 2e-6)
all_feature$noncoding_mutation_rate[index_largenoncoding_mutation_rate]<-(all_feature$noncoding_mutation_rate[index_largenoncoding_mutation_rate]- 4e-6)/10 + 4e-6
all_feature$noncoding_mutation_rate[index_lownoncoding_mutation_rate]<-(all_feature$noncoding_mutation_rate[index_lownoncoding_mutation_rate]- 2e-6)/10 + 2e-6

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
index_pOG<-which(prediction$OG_probability>OG_threshold & prediction$TSG_probability<TSG_threshold & OG_CGC!="1")
index_pTSG<-which(prediction$TSG_probability>TSG_threshold & prediction$OG_probability<OG_threshold & TSG_CGC!="1")
index_TSG<-which((TSG_CGC=="1" & OG_CGC!="1"))
index_OG<-which((OG_CGC=="1" & TSG_CGC!="1"))

index<-rep("",nrow(anno))
index[which(NG=="1")]<-"NG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

all_feature<-cbind(all_feature,index)
allgene2<-all_feature[all_feature$index!="" & !is.na(all_feature$noncoding_mutation_rate),]
allgene2$index<-factor(allgene2$index,levels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))
pdf("Raw_figures/Figure_S4F_noncoding_mutation_evaluation.pdf", family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-TSG","CGC-OG"),c("Novel DORGE-TSG","Novel DORGE-OG"),c("CGC-TSG","NG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("Novel DORGE-TSG","NG"),c("Novel DORGE-TSG","CGC-OG"),c("CGC-TSG","Novel DORGE-OG"));
ggboxplot(data = allgene2,x = "index", y="noncoding_mutation_rate",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_x_discrete(labels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))+scale_y_continuous(breaks = c(2e-6,4e-6,6e-6),labels = c("2e-6", "4e-6", "2.2e-5"))+ labs(x = "", y = "Noncoding mutation rate") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "l"),comparisons = my_comparisons,size=2.5)
garbage <-dev.off()

######### Figure S4C: shRNA screening data evaluation for DORGE prediction #########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_4");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("data/shRNA_cell_proliferation_rates.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","Cell_proliferation_rate")]
all_feature$Cell_proliferation_rate <- as.numeric(all_feature$Cell_proliferation_rate)
index_low_signal<-which(all_feature$Cell_proliferation_rate< -2)
all_feature$Cell_proliferation_rate[index_low_signal]<-(all_feature$Cell_proliferation_rate[index_low_signal]+2)/10-2

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")

index_pOG<-which(prediction$OG_probability> OG_threshold & prediction$TSG_probability<TSG_threshold & OG_CGC!="1")
index_pTSG<-which(prediction$TSG_probability>TSG_threshold & prediction$OG_probability< OG_threshold & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")

index<-rep("",nrow(anno))
index[which(NG=="1")]<-"NG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"
all_feature<-cbind(all_feature,index)
allgene2<-all_feature[all_feature$index!="" & !is.na(all_feature$Cell_proliferation_rate),]
allgene2$index<-factor(allgene2$index,levels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))

pdf("Raw_figures/Figure_S4C_shRNA_evaluation.pdf", family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-OG","CGC-TSG"),c("Novel DORGE-OG","Novel DORGE-TSG"),c("CGC-TSG","NG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("Novel DORGE-TSG","NG"),c("Novel DORGE-OG","CGC-TSG"),c("CGC-OG","Novel DORGE-TSG"));
ggboxplot(data = allgene2,x = "index", y="Cell_proliferation_rate",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_x_discrete(labels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))+scale_y_continuous(breaks = c(-2,-1,0,1,2),labels = c("-2", "-1", "0", "1", "2"))+ labs(x = "", y = "Cell proliferation rate") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "l"),comparisons = my_comparisons,size=2.5)
garbage <-dev.off()

######## Figure S4D Exoncons evaluation########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_4");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","Exon_conservation_phastCons_score")]
all_feature$Exon_conservation_phastCons_score <- as.numeric(all_feature$Exon_conservation_phastCons_score)
index_largeExon_conservation_phastCons_score<-which(all_feature$Exon_conservation_phastCons_score<0.4)
all_feature$Exon_conservation_phastCons_score[index_largeExon_conservation_phastCons_score]<-(all_feature$Exon_conservation_phastCons_score[index_largeExon_conservation_phastCons_score]- 0.4)/2 + 0.4
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")

index_pOG<-which(prediction$OG_probability> OG_threshold & prediction$TSG_probability<TSG_threshold & OG_CGC!="1")
index_pTSG<-which(prediction$TSG_probability>TSG_threshold & prediction$OG_probability< OG_threshold & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")

index<-rep("",nrow(anno))
index[which(NG=="1")]<-"NG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

all_feature<-cbind(all_feature,index)
allgene2<-all_feature[all_feature$index!="",]
allgene2$index<-factor(allgene2$index,levels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))
pdf("Raw_figures/Figure_S4D_phastcons_evaluation.pdf", family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-TSG","CGC-OG"),c("Novel DORGE-TSG","Novel DORGE-OG"),c("CGC-TSG","NG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("Novel DORGE-TSG","NG"),c("Novel DORGE-TSG","CGC-OG"),c("CGC-TSG","Novel DORGE-OG"));

ggboxplot(data = allgene2,x = "index", y="Exon_conservation_phastCons_score",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_x_discrete(labels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))+scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1),labels = c("0.3", "0.4", "0.6", "0.8","1.0"))+ labs(x = "", y = "phastCons") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <-dev.off()

######## Figure 4E phyloP evaluation########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_4");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("data/phyloP46_processed.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","exon_phyloP")]
all_feature$exon_phyloP <- as.numeric(all_feature$exon_phyloP)
index_largeexon_phyloP<-which(all_feature$exon_phyloP>3)
index_lowexon_phyloP<-which(all_feature$exon_phyloP<0)
all_feature$exon_phyloP[index_largeexon_phyloP]<-(all_feature$exon_phyloP[index_largeexon_phyloP]-3)/10+3
all_feature$exon_phyloP[index_lowexon_phyloP]<-(all_feature$exon_phyloP[index_lowexon_phyloP])/10


TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
index_pOG<-which(prediction$OG_probability> OG_threshold & prediction$TSG_probability<TSG_threshold & OG_CGC!="1")
index_pTSG<-which(prediction$TSG_probability>TSG_threshold & prediction$OG_probability< OG_threshold & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")

index<-rep("",nrow(anno))
index[which(NG=="1")]<-"NG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

all_feature<-cbind(all_feature,index)
allgene2<-all_feature[all_feature$index!="",]
allgene2$index<-factor(allgene2$index,levels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))
pdf("Raw_figures/Figure_4E_PhyloP_evaluation.pdf", family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-TSG","CGC-OG"),c("Novel DORGE-TSG","Novel DORGE-OG"),c("CGC-TSG","NG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("Novel DORGE-TSG","NG"),c("Novel DORGE-TSG","CGC-OG"),c("CGC-TSG","Novel DORGE-OG"));

ggboxplot(data = allgene2,x = "index", y="exon_phyloP",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_x_discrete(labels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))+scale_y_continuous(breaks = c(-1,0,1,2,3,4),labels = c("-1","0", "1", "2", "3","4"))+ labs(x = "", y = "phyloP score") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <-dev.off()


######## Figure S4E ncRVIS evaluation ########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_4");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","ncRVIS_score")]
all_feature$ncRVIS_score <- as.numeric(all_feature$ncRVIS_score)
index_large_ncRVIS_score<-which(all_feature$ncRVIS_score > 2)
all_feature$ncRVIS_score[index_large_ncRVIS_score]<-(all_feature$ncRVIS_score[index_large_ncRVIS_score]-2)/20+2
index_low_ncRVIS_score<-which(all_feature$ncRVIS_score < -2)
all_feature$ncRVIS_score[index_low_ncRVIS_score]<-(all_feature$ncRVIS_score[index_low_ncRVIS_score]+2)/10-2

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
index_pOG<-which(prediction$OG_probability> OG_threshold & prediction$TSG_probability<TSG_threshold & OG_CGC!="1")
index_pTSG<-which(prediction$TSG_probability>TSG_threshold & prediction$OG_probability< OG_threshold & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")

index<-rep("",nrow(anno))
index[which(NG=="1")]<-"NG"
index[index_pOG]<-"Novel DORGE-OG"
index[index_pTSG]<-"Novel DORGE-TSG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"
all_feature<-cbind(all_feature,index)
allgene2<-all_feature[all_feature$index!="",]
allgene2$index<-factor(allgene2$index,levels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))
pdf("Raw_figures/Figure_S4E_ncRVIS_evaluation.pdf", family="ArialMT", width=2.275, height=3)
my_comparisons <- list(c("CGC-OG","CGC-TSG"),c("Novel DORGE-OG","Novel DORGE-TSG"),c("CGC-TSG","NG"),c("Novel DORGE-OG","NG"),c("CGC-OG","NG"),c("Novel DORGE-TSG","NG"),c("Novel DORGE-OG","CGC-TSG"),c("CGC-OG","Novel DORGE-TSG"));

ggboxplot(data = allgene2,x = "index", y="ncRVIS_score",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+scale_x_discrete(labels=c("CGC-OG", "CGC-TSG", "NG", "Novel DORGE-OG", "Novel DORGE-TSG"))+scale_y_continuous(breaks = c(-2,0,2),labels = c("-2", "0", "2"))+ labs(x = "", y = "ncRVIS") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "l"),comparisons = my_comparisons,size=2.5)
garbage <-dev.off()