####################
##### Boxplots #####
####################

###### The boxplots for different gene categories to show the distribution of top-ranking features used in DORGE-TSG prediction.
###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/genebody_canyon_genes.txt: Compiled gene-body canyon genes from Canyon_with_genebody_processed.bed at Epigenetics_processing folder from the DORGE_pipeline project https://github.com/biocq/DORGE_pipeline
###### Output: Figure_1C_H3K4me3_peak_length.pdf: Box plots showing the distribution of tri-methylation on histone H3 lysine 4 (H3K4me3) mean peak length among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_1D_VEST_score.pdf: Box plots showing the distribution of Variant Effect Scoring Tool (VEST) score among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_1E_Missense_damaging_benign_ratio.pdf: Box plots showing the distribution of missense damaging/benign ratio among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_1F_Missense_entropy.pdf: Box plots showing the distribution of missense entropy among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_1G_pLI_score.pdf: Box plots showing the distribution of pLI score among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_1H_super_enhancer_percentage.pdf: Box plots showing the distribution of super enhancer percentage among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_1I_gene-body_hypermethylation.pdf: Box plots showing the distribution of gene-body hypermethylation among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1A_Height_of_H3K4me3_peaks.pdf: Box plots showing the distribution of tri-methylation of histone H3 lysine 4 (H3K4me3) peak height among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1B_H3K79me2_peak_length.pdf: Box plots showing the distribution of di-methylation of histone H3 lysine 79 (H3K79me2) peak length among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1C_Height_of_H3K79me2_peaks.pdf: Box plots showing the distribution of di-methylation of histone H3 lysine 79 (H3K79me2) peak height among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1D_NonSilent_silent_ratio.pdf: Box plots showing the distribution of non-silent/silent ratio among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1E_LoF_o_e_constraint.pdf: Box plots showing the distribution of LoF o/e constraint among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1F_Exon_conservation_phastCons_score.pdf: Box plots showing the distribution of phastCons score among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1G_ncGERP_score.pdf: Box plots showing the distribution of non-coding Genomic Evolutionary Rate Profiling (ncGERP) score among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1H_PolyPhen_2_score.pdf: Box plots showing the distribution of PolyPhen-2 score among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1I_H3K36me3_peak_length.pdf: Box plots showing the distribution of tri-methylation of histone H3 lysine 36 (H3K36me3) peak length among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1J_H4K20me1_peak_length.pdf: Box plots showing the distribution of mono-methylation of histone H4 lysine 20 (H4K20me1) peak length among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1K_H3K9ac_peak_length.pdf: Box plots showing the distribution of histone H3 lysine 9 acetylation (H3K9ac) peak length among CGC-OG, CGC-TSG, and neutral gene (NG) sets
###### Output: Figure_S1M_gene-body_hypermethylation_at_canyon_genes.pdf: Box plots showing the distribution of genebody differential methylation at genebody canyon genes


options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggpubr","scales")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

suppressMessages(library("ggpubr"))
suppressMessages(library("scales"))

################## H3K4me3 peak length ################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","H3K4me3_peak_length")]
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>4000)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-4000)/4+4000
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))
pdf("Raw_figures/Figure_1C_H3K4me3_peak_length.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="H3K4me3_peak_length",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,2000,4000,6000),labels = c("0", "2000", "4000", "12000"))+ labs(x = "", y = "H3K4me3 peak length") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################## H3K4me3 peak height ################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Height_of_H3K4me3_peaks")]
allgene$Height_of_H3K4me3_peaks <- as.numeric(allgene$Height_of_H3K4me3_peaks)
index_higherK4me3<-which(allgene$Height_of_H3K4me3_peaks>1000)
allgene$Height_of_H3K4me3_peaks[index_higherK4me3]<-(allgene$Height_of_H3K4me3_peaks[index_higherK4me3]-1000)/10+1000
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))
pdf("Raw_figures/Figure_S1A_Height_of_H3K4me3_peaks.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="Height_of_H3K4me3_peaks",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,500,1000,1200),labels = c("0", "500", "1000","3000"))+ labs(x = "", y = "H3K4me3 peak height") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################# H3K79me2 peak length #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","H3K79me2_peak_length")]
allgene$H3K79me2_peak_length <- as.numeric(allgene$H3K79me2_peak_length)
index_BroadK79me2<-which(allgene$H3K79me2_peak_length>4000)
allgene$H3K79me2_peak_length[index_BroadK79me2]<-(allgene$H3K79me2_peak_length[index_BroadK79me2]-4000)/4+4000
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))
pdf("Raw_figures/Figure_S1B_H3K79me2_peak_length.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="H3K79me2_peak_length",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,2000,4000,6000),labels = c("0", "2000", "4000", "12000"))+ labs(x = "", y = "H3K79me2 peak length") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################## H3K79me2 peak height ################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Height_of_H3K79me2_peaks")]
allgene$Height_of_H3K79me2_peaks <- as.numeric(allgene$Height_of_H3K79me2_peaks)
index_higherH3K79me2<-which(allgene$Height_of_H3K79me2_peaks>40)
allgene$Height_of_H3K79me2_peaks[index_higherH3K79me2]<-(allgene$Height_of_H3K79me2_peaks[index_higherH3K79me2]-40)/10+40
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))
pdf("Raw_figures/Figure_S1C_Height_of_H3K79me2_peaks.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="Height_of_H3K79me2_peaks",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,20,40,50),labels = c("0", "20", "40","140"))+ labs(x = "", y = "H3K79me2 peak height") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################# H3K36me3 peak length #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","H3K36me3_peak_length")]
allgene$H3K36me3_peak_length <- as.numeric(allgene$H3K36me3_peak_length)
index_BroadH3K36me3<-which(allgene$H3K36me3_peak_length>4000)
allgene$H3K36me3_peak_length[index_BroadH3K36me3]<-(allgene$H3K36me3_peak_length[index_BroadH3K36me3]-4000)/4+4000
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))
pdf("Raw_figures/Figure_S1I_H3K36me3_peak_length.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="H3K36me3_peak_length",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,2000,4000,6000),labels = c("0", "2000", "4000", "12000"))+ labs(x = "", y = "H3K36me3 peak length") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# H4K20me1 peak length #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","H4K20me1_peak_length")]
allgene$H4K20me1_peak_length <- as.numeric(allgene$H4K20me1_peak_length)
index_BroadH4K20me1<-which(allgene$H4K20me1_peak_length>4000)
allgene$H4K20me1_peak_length[index_BroadH4K20me1]<-(allgene$H4K20me1_peak_length[index_BroadH4K20me1]-4000)/4+4000
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))
pdf("Raw_figures/Figure_S1J_H4K20me1_peak_length.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="H4K20me1_peak_length",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,2000,4000,6000),labels = c("0", "2000", "4000", "12000"))+ labs(x = "", y = "H4K20me1 peak length") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# H3K9ac peak length #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","H3K9ac_peak_length")]
allgene$H3K9ac_peak_length <- as.numeric(allgene$H3K9ac_peak_length)
index_BroadH3K9ac<-which(allgene$H3K9ac_peak_length>4000)
allgene$H3K9ac_peak_length[index_BroadH3K9ac]<-(allgene$H3K9ac_peak_length[index_BroadH3K9ac]-4000)/4+4000
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))
pdf("Raw_figures/Figure_S1K_H3K9ac_peak_length.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="H3K9ac_peak_length",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,2000,4000,6000),labels = c("0", "2000", "4000", "12000"))+ labs(x = "", y = "H3K9ac peak length") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# VEST score #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","VEST_score")]


TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index<-rep("",nrow(allgene))
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index[which(NG=="1")]<-"NG"
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_1D_VEST_score.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-TSG","CGC-OG"))
ggboxplot(data = allgene2,x = "index", y="VEST_score",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "VEST score")+ theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# Exon conservation phastCons score #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Exon_conservation_phastCons_score")]

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index[index_OG]<-"CGC-OG"
index[index_TSG]<-"CGC-TSG"
allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_S1F_Exon_conservation_phastCons_score.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-TSG","CGC-OG"))
ggboxplot(data = allgene2,x = "index", y="Exon_conservation_phastCons_score",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "phastCons score")+ theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################# Gene-body differential hypermethylation at gene-body canyon genes #################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Gene_body_hypermethylation_in_cancer","Gene_body_canyon_hypermethylation_in_cancer")]

genebody_canyon_genes <- read.table("data/genebody_canyon_genes.txt", header=F)
canyon_methyl_diff <- read.table("data/Canyon_methylation.txt", header=T, sep="\t",fill=TRUE,quote = "")
genebody_canyon_hypermethylation_diff<-as.numeric(as.character(canyon_methyl_diff$Cancer_median))-as.numeric(as.character(canyon_methyl_diff$Normal_median))
allgene$genebody_canyon_hypermethylation_diff<-genebody_canyon_hypermethylation_diff

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene<-allgene[which(allgene$Gene%in%genebody_canyon_genes$V1),]
allgene2<-allgene[allgene$index!="" & allgene$Gene%in%genebody_canyon_genes$V1,]

my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-OG","CGC-TSG"))

#ggboxplot(data = allgene2,x = "index", y="Gene_body_canyon_hypermethylation_in_cancer",color="black",outlier.size = 0.6,width=0.7, fill="index",palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "Gene-body hypermethylation ") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
pdf("Raw_figures/Figure_S1M_gene-body_hypermethylation_at_canyon_genes.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG", "NG"),c("NG","CGC-TSG"),c("CGC-OG", "CGC-TSG"))
ggboxplot(data = allgene2,x = "index", y="genebody_canyon_hypermethylation_diff",color="black",outlier.size = 0.6,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "Gene-body differential methylation\nat gene-body canyon genes") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# Gene-body differential hypermethylation #################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Gene_body_hypermethylation_in_cancer")]

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-OG","CGC-TSG"))

pdf("Raw_figures/Figure_1I_gene-body_hypermethylation.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG", "NG"),c("NG","CGC-TSG"),c("CGC-OG", "CGC-TSG"))
ggboxplot(data = allgene2,x = "index", y="Gene_body_hypermethylation_in_cancer",color="black",outlier.size = 0.6,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "Gene-body hypermethylation") + theme_bw() + theme(axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# Missense entropy #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Missense_entropy")]

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

# John and Draper's modulus transformation
modulus_trans <- function(lambda){
   trans_new("modulus",
     transform = function(y){
        if(lambda != 0){
           yt <- sign(y) * (((abs(y) + 1) ^ lambda - 1) / lambda)
        }
        return(yt)
     },
     inverse = function(yt){
        if(lambda != 0){
           y <- ((abs(yt) * lambda + 1)  ^ (1 / lambda) - 1) * sign(yt)
        } else {
           y <- (exp(abs(yt)) - 1) * sign(yt)
           
        }
        return(y)
     }
   )
}
pdf("Raw_figures/Figure_1F_Missense_entropy.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-OG","CGC-TSG"))
ggboxplot(data = allgene2,x = "index", y="Missense_entropy",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+labs(x = "", y = "Missense entropy") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)+ guides(fill=FALSE)+ scale_y_continuous(trans = modulus_trans(-0.7),breaks = c(0,0.2,0.6,1,4))
garbage <- dev.off()


################# Missense damaging/benign ratio #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Missense_damaging_benign_ratio")]

allgene$Missense_damaging_benign_ratio <- as.numeric(allgene$Missense_damaging_benign_ratio)
index_higher<-which(allgene$Missense_damaging_benign_ratio>2)
allgene$Missense_damaging_benign_ratio[index_higher]<-(allgene$Missense_damaging_benign_ratio[index_higher]-2)/50+2

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_1E_Missense_damaging_benign_ratio.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-OG","CGC-TSG"))
ggboxplot(data = allgene2,x = "index", y="Missense_damaging_benign_ratio",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,1,2,3,4),labels = c("0", "1", "2","52","102"))+ labs(x = "", y = "Missense damaging/benign ratio") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################# Super_enhancer_percentage #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","Super_enhancer_percentage")]
allgene$Super_enhancer_percentage <- as.numeric(allgene$Super_enhancer_percentage)
index_higher<-which(allgene$Super_enhancer_percentage>0.2)
allgene$Super_enhancer_percentage[index_higher]<-(allgene$Super_enhancer_percentage[index_higher]-0.2)/8+0.2

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_1H_super_enhancer_percentage.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-TSG","CGC-OG"))
ggboxplot(data = allgene2,x = "index", y="Super_enhancer_percentage",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+  scale_y_continuous(breaks = c(0,0.05,0.10,0.15,0.20,0.25),labels = c("0","0.05", "0.10", "0.15","0.20","0.60"))+ labs(x = "", y = "Super enhancer percentage") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# ncGERP score #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","ncGERP_score")]

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_S1G_ncGERP_score.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-OG","CGC-TSG"))
ggboxplot(data = allgene2,x = "index", y="ncGERP_score",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "ncGERP score") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()

################# PolyPhen-2 score #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","PolyPhen_2_score")]

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_S1H_PolyPhen_2_score.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG","NG"),c("CGC-TSG","NG"),c("CGC-OG","CGC-TSG"))
ggboxplot(data = allgene2,x = "index", y="PolyPhen_2_score",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "PolyPhen-2 score") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################# pLI score #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","pLI_score")]

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_1G_pLI_score.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-OG", "CGC-TSG"))
ggboxplot(data = allgene2,x = "index", y="pLI_score",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "pLI score") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)
garbage <- dev.off()


################# LoF o/e constraint #################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","LoF_o_e_constraint")]
allgene$LoF_o_e_constraint <- as.numeric(allgene$LoF_o_e_constraint)
index_higher<-which(allgene$LoF_o_e_constraint>8)
allgene$LoF_o_e_constraint[index_higher]<-(allgene$LoF_o_e_constraint[index_higher]-8)/10+8
index_lower<-which(allgene$LoF_o_e_constraint< -2)
allgene$LoF_o_e_constraint[index_lower]<-(allgene$LoF_o_e_constraint[index_lower]+2)/10-2

TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_S1E_LoF_o_e_constraint.pdf", family="ArialMT", width=1.6, height=3)
my_comparisons <- list(c("CGC-OG", "NG"),c("CGC-TSG","NG"),c("CGC-TSG", "CGC-OG"))

ggboxplot(data = allgene2,x = "index", y="LoF_o_e_constraint",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "LoF o/e constraint") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)+  scale_y_continuous(breaks = c(-2,0,2,4,8),labels = c("-2", "0", "2", "4", "8"))
garbage <- dev.off()

################# NonSilent/silent ratio #################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c("Gene","NonSilent_silent_ratio")]
allgene$NonSilent_silent_ratio <- as.numeric(allgene$NonSilent_silent_ratio)
index_higher<-which(allgene$NonSilent_silent_ratio>5)
allgene$NonSilent_silent_ratio[index_higher]<-(allgene$NonSilent_silent_ratio[index_higher]-5)/100+5
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")

index<-rep("",nrow(allgene))
index[which(NG=="1")]<-"NG"
index[index_TSG]<-"CGC-TSG"
index[index_OG]<-"CGC-OG"

allgene<-cbind(allgene,index)
allgene2<-allgene[allgene$index!="",]

pdf("Raw_figures/Figure_S1D_NonSilent_silent_ratio.pdf", family="ArialMT", width=1.6, height=3)
ggboxplot(data = allgene2,x = "index", y="NonSilent_silent_ratio",color="black",outlier.size = 0.2,width=0.7, fill="index",lwd=0.25,palette = c("#E41A1C","#377EB8","#CCCCCC"))+ labs(x = "", y = "Non-silent/silent ratio") + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)+ scale_y_continuous(breaks = c(0,1,3,5,5.9,7.95),labels = c("0", "1", "3", "5", "100","300"))
garbage <- dev.off()
