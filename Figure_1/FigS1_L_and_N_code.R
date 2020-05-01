#####################
##### Barplots  #####
#####################

###### The barplots showing the enrichment (-log10(P-values of Fisher's exact test)) of Broad H3K4me3 genes and gene-body hypermethylation canyon genes in different gene categories.
###### Input: ../Gene_set_new.txt: Gene annotation file
###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: data/genebody_canyon_genes.txt: Gene-body canyon genes (Compiled from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1492-3/figures/5)
###### Output: Figure_S1L_Broad_K4me3_enrichment.pdf: Fishers exact test for genes with broad H3K4me3 peaks (length > 4,000 bp) in the CGC-OG, CGC-TSG, and NG sets.
###### Output: Figure_S1N_genebody_canyon_enrichment.pdf: Fishers exact test for genes with hypermethylated gene-body canyon in the CGC-OG, CGC-TSG, and NG sets.
options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggplot2")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("ggplot2"))

################################### Broad H3K4me3 peak gene enrichment ###################################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c(1,66)]
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>4000)
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-4000)/4+4000
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_TSG<-which(TSG_CGC=="1" & OG_CGC!="l")
index_OG<-which(OG_CGC=="1" & TSG_CGC!="1")
index_NG<-which(NG=="1")

known_tsg_genes<-as.character(allgene[index_TSG,1])
known_og_genes<-as.character(allgene[index_OG,1])
neutral_genes<-as.character(allgene[index_NG,1])

nonknownTSG<-setdiff(allgene[,1],known_tsg_genes)
nonknownOG<-setdiff(allgene[,1],known_og_genes)
nonNG<-setdiff(allgene[,1],neutral_genes)

Broad_K4me3_genes<-as.character(allgene$Gene[allgene$H3K4me3_peak_length>4000])
nonBroad_K4me3_genes<-setdiff(allgene[,1],Broad_K4me3_genes)

norm_number<-min(c(length(known_tsg_genes),length(known_og_genes),length(neutral_genes)))

#CGC-TSG
CGC_TSG_Broad_K4me3<-as.character(known_tsg_genes[which(known_tsg_genes%in%Broad_K4me3_genes)])
CGC_TSG_nonBroad_K4me3<-as.character(known_tsg_genes[which(known_tsg_genes%in%nonBroad_K4me3_genes)])
nonCGC_TSG_Broad_K4me3<-as.character(nonknownTSG[which(nonknownTSG%in%Broad_K4me3_genes)])
nonCGC_TSG_nonBroad_K4me3<-as.character(nonknownTSG[which(nonknownTSG%in%nonBroad_K4me3_genes)])
fisher<-matrix(c(length(CGC_TSG_Broad_K4me3),length(CGC_TSG_nonBroad_K4me3),length(nonCGC_TSG_Broad_K4me3), length(nonCGC_TSG_nonBroad_K4me3))/(length(known_tsg_genes)/norm_number),nrow = 2,dimnames =list(c("Broad_K4me3_genes", "nonBroad_K4me3_genes"),c("TSG", "nonTSG")))
pval_CGC_TSG_Broad_K4me3<-fisher.test(fisher, alternative = "g")$p.value


CGC_OG_Broad_K4me3<-as.character(known_og_genes[which(known_og_genes%in%Broad_K4me3_genes)])
CGC_OG_nonBroad_K4me3<-as.character(known_og_genes[which(known_og_genes%in%nonBroad_K4me3_genes)])
nonCGC_OG_Broad_K4me3<-as.character(nonknownOG[which(nonknownOG%in%Broad_K4me3_genes)])
nonCGC_OG_nonBroad_K4me3<-as.character(nonknownOG[which(nonknownOG%in%nonBroad_K4me3_genes)])
fisher<-matrix(c(length(CGC_OG_Broad_K4me3),length(CGC_OG_nonBroad_K4me3),length(nonCGC_OG_Broad_K4me3), length(nonCGC_OG_nonBroad_K4me3))/(length(known_og_genes)/norm_number),nrow = 2,dimnames =list(c("Broad_K4me3_genes", "nonBroad_K4me3_genes"),c("OG", "nonOG")))
pval_CGC_OG_Broad_K4me3<-fisher.test(fisher, alternative = "g")$p.value


NG_Broad_K4me3<-as.character(neutral_genes[which(neutral_genes%in%Broad_K4me3_genes)])
NG_nonBroad_K4me3<-as.character(neutral_genes[which(neutral_genes%in%nonBroad_K4me3_genes)])
nonNG_Broad_K4me3<-as.character(nonNG[which(nonNG%in%Broad_K4me3_genes)])
nonNG_nonBroad_K4me3<-as.character(nonNG[which(nonNG%in%nonBroad_K4me3_genes)])
fisher<-matrix(c(length(NG_Broad_K4me3),length(NG_nonBroad_K4me3),length(nonNG_Broad_K4me3), length(nonNG_nonBroad_K4me3))/(length(neutral_genes)/norm_number),nrow = 2,dimnames =list(c("Broad_K4me3_genes", "nonBroad_K4me3_genes"),c("NG", "nonNG")))
pval_NG_Broad_K4me3<-fisher.test(fisher, alternative = "g")$p.value

pval<-c(pval_CGC_TSG_Broad_K4me3,pval_CGC_OG_Broad_K4me3,pval_NG_Broad_K4me3-0.1)#add Pseudo-count
pval<- -log10(pval)
Type<-c("CGC-TSG","CGC-OG","NG")
dat<-data.frame(Type,pval)
dat$Type<-factor(dat$Type,levels=c("CGC-TSG","CGC-OG","NG"))

pdf("Raw_figures/Figure_S1L_Broad_K4me3_enrichment.pdf", family="ArialMT", width=1.6, height=3)
ggplot(data=dat, aes(x=as.factor(Type), y=pval,fill=as.factor(Type)))+geom_bar(stat="identity", alpha=1, position="dodge",width=0.85) +labs(x="",y=expression("-log"[10]~"Enrichment"~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.title.x = element_text(hjust=1),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(colour = "black"),legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ scale_fill_manual(values=rev(c("#CCCCCC","#E41A1C","#377EB8")))
garbage <- dev.off()



################################### Gene-body canyon gene enrichment ###################################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_1");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-all_feature[,c(1,66)]

genebody_canyon_genes <- read.table("data/genebody_canyon_genes.txt", header=F)
genebody_canyon_genes <- genebody_canyon_genes$V1
non_genebody_canyon_genes<-setdiff(allgene[,1],genebody_canyon_genes)
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

index_TSG<-which(TSG_CGC=="1" & OG_CGC!="1")
index_OG<-which(OG_CGC=="1" & TSG_CGC!="l")
index_NG<-which(NG=="1")

known_tsg_genes<-as.character(allgene[index_TSG,1])
known_og_genes<-as.character(allgene[index_OG,1])
neutral_genes<-as.character(allgene[index_NG,1])

nonknownTSG<-setdiff(allgene[,1],known_tsg_genes)
nonknownOG<-setdiff(allgene[,1],known_og_genes)
nonNG<-setdiff(allgene[,1],neutral_genes)


norm_number<-min(c(length(known_tsg_genes),length(known_og_genes),length(neutral_genes)))

#CGC-TSG
CGC_TSG_genebody_canyon<-as.character(known_tsg_genes[which(known_tsg_genes%in%genebody_canyon_genes)])
CGC_TSG_non_genebody_canyon<-as.character(known_tsg_genes[which(known_tsg_genes%in%non_genebody_canyon_genes)])
nonCGC_TSG_genebody_canyon<-as.character(nonknownTSG[which(nonknownTSG%in%genebody_canyon_genes)])
nonCGC_TSG_non_genebody_canyon<-as.character(nonknownTSG[which(nonknownTSG%in%non_genebody_canyon_genes)])
fisher<-matrix(c(length(CGC_TSG_genebody_canyon),length(CGC_TSG_non_genebody_canyon),length(nonCGC_TSG_genebody_canyon), length(nonCGC_TSG_non_genebody_canyon))/(length(known_tsg_genes)/norm_number),nrow = 2,dimnames =list(c("genebody_canyon_genes", "non_genebody_canyon_genes"),c("TSG", "nonTSG")))
pval_CGC_TSG_genebody_canyon<-fisher.test(fisher, alternative = "g")$p.value


CGC_OG_genebody_canyon<-as.character(known_og_genes[which(known_og_genes%in%genebody_canyon_genes)])
CGC_OG_non_genebody_canyon<-as.character(known_og_genes[which(known_og_genes%in%non_genebody_canyon_genes)])
nonCGC_OG_genebody_canyon<-as.character(nonknownOG[which(nonknownOG%in%genebody_canyon_genes)])
nonCGC_OG_non_genebody_canyon<-as.character(nonknownOG[which(nonknownOG%in%non_genebody_canyon_genes)])
fisher<-matrix(c(length(CGC_OG_genebody_canyon),length(CGC_OG_non_genebody_canyon),length(nonCGC_OG_genebody_canyon), length(nonCGC_OG_non_genebody_canyon))/(length(known_og_genes)/norm_number),nrow = 2,dimnames =list(c("genebody_canyon_genes", "non_genebody_canyon_genes"),c("OG", "nonOG")))
pval_CGC_OG_genebody_canyon<-fisher.test(fisher, alternative = "g")$p.value


NG_genebody_canyon<-as.character(neutral_genes[which(neutral_genes%in%genebody_canyon_genes)])
NG_non_genebody_canyon<-as.character(neutral_genes[which(neutral_genes%in%non_genebody_canyon_genes)])
nonNG_genebody_canyon<-as.character(nonNG[which(nonNG%in%genebody_canyon_genes)])
nonNG_non_genebody_canyon<-as.character(nonNG[which(nonNG%in%non_genebody_canyon_genes)])
fisher<-matrix(c(length(NG_genebody_canyon),length(NG_non_genebody_canyon),length(nonNG_genebody_canyon), length(nonNG_non_genebody_canyon))/(length(neutral_genes)/norm_number),nrow = 2,dimnames =list(c("genebody_canyon_genes", "non_genebody_canyon_genes"),c("NG", "nonNG")))
pval_NG_genebody_canyon<-fisher.test(fisher, alternative = "g")$p.value

pval<-c(pval_CGC_TSG_genebody_canyon,pval_CGC_OG_genebody_canyon,pval_NG_genebody_canyon-0.1)#add Pseudo-count
pval<- -log10(pval)
Type<-c("CGC-TSG","CGC-OG","NG")
dat<-data.frame(Type,pval)
dat$Type<-factor(dat$Type,levels=c("CGC-TSG","CGC-OG","NG"))

pdf("Raw_figures/Figure_S1N_genebody_canyon_enrichment.pdf", family="ArialMT", width=1.6, height=3)
ggplot(data=dat, aes(x=as.factor(Type), y=pval,fill=as.factor(Type)))+geom_bar(stat="identity", alpha=1, position="dodge",width=0.85)+labs(x="",y=expression("-log"[10]~"Enrichment"~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.title.x = element_text(hjust=1),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(colour = "black"),legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ scale_fill_manual(values=rev(c("#CCCCCC","#E41A1C","#377EB8")))
garbage <- dev.off()