####################
##### Barplots #####
####################

###### The enrichment of Protein-Protein interaction (PPI) network hub genes in DORGE-predicted TSG/OGs, NGs, CGC genes and other functional genomics gene sets (or features).

###### Input: ../Gene_set_new.txt: Gene annotation file
###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/ER_genes.txt: Epigenetic regulator (ER) gene set (Compiled from https://genome.cshlp.org/content/29/4/532.long and https://academic.oup.com/database/article/doi/10.1093/database/bav067/2433200)
###### Input: data/gene_essentiality/ essential_genes.txt and non_essential_genes.txt: Essential and non-essential gene lists (from OGEE database)
###### Input: data/genebody_canyon_genes.txt: Gene-body canyon genes (Compiled from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1492-3)
###### Input: data/genebody_differential_methyl_COSMIC.txt (Gene-body methylation data in cancer and normal samples compiled from the v90 methylation data in COSMIC website)
###### Input: data/HKG_processed.txt: Housekeeping gene (HKG) list that was downloaded from https://www.tau.ac.il/~elieis/HKG/

###### Output: Figure_5C_BioGRID_network_hub_enrichment.pdf: PPI network hub gene enrichment in DORGE-predicted TSG/OGs, NGs and CGC genes
###### Output: Figure_5D_BioGRID_network_hub_enrichment.pdf: PPI network hub gene enrichment in functional genomics gene sets (or features)


options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("plyr","ggsignif","ggpubr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

suppressMessages(library("ggpubr"))
suppressMessages(library(ggsignif))

TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005
######### Figure 5C and 5D: Gene enrichment in PPI network #########

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_5");
anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
anno<-cbind(anno[,c(1,4,5,6)],prediction[,2:3])
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","BioGRID_betweenness","BioGRID_clossness","Log_BioGRID_degree","H3K4me3_peak_length","LoF_mutations_kb","Missense_mutations_kb","Missense_o_e_constraint","LoF_o_e_constraint")]
all_feature$BioGRID_betweenness <- as.numeric(all_feature$BioGRID_betweenness)
all_feature$BioGRID_clossness <- as.numeric(all_feature$BioGRID_clossness)
all_feature$Log_BioGRID_degree <- as.numeric(all_feature$Log_BioGRID_degree)
all_feature$H3K4me3_peak_length <- as.numeric(all_feature$H3K4me3_peak_length)
all_feature$LoF_mutations_kb <- as.numeric(all_feature$LoF_mutations_kb)
all_feature$Missense_o_e_constraint <- as.numeric(all_feature$Missense_o_e_constraint)
all_feature$LoF_o_e_constraint <- as.numeric(all_feature$LoF_o_e_constraint)

index_high<-which(all_feature$BioGRID_betweenness>1e-3)
all_feature$BioGRID_betweenness[index_high]<-(all_feature$BioGRID_betweenness[index_high]-1e-3)/10+1e-3

index_high<-which(all_feature$BioGRID_clossness>0.4)
all_feature$BioGRID_clossness[index_high]<-(all_feature$BioGRID_clossness[index_high]-0.4)/10+0.4

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

quantile_table<-quantile(dat2$Log_BioGRID_degree, prob = seq(0, 1, length = 21), type = 5,na.rm=T) #95% percentile
index_high<-which(dat2$Log_BioGRID_degree>as.numeric(quantile_table[length(quantile_table)-1]))
hub_genes<-as.character(dat2[index_high,1])
nonhub_genes<-as.character(dat2[-index_high,1])

essential <- read.table("data/gene_essentiality/essential_genes.txt", header=F, sep="\t")
essential<-as.character(essential[,1])
non_essential <- read.table("data/gene_essentiality/non_essential_genes.txt", header=F, sep="\t")
non_essential<-as.character(non_essential[,1])

ER <- read.table("data/ER_genes.txt", header=F, sep="\t")
ER_genes<-as.character(ER[,1])
non_ER_genes<-as.character(anno$Gene[-which(anno$Gene%in%ER_genes)])

canyon_genes <- read.table("data/genebody_canyon_genes.txt", header=F)
canyon_genes <- as.character(canyon_genes[,1])
noncanyon_genes<-as.character(anno$Gene[-which(anno$Gene%in%canyon_genes)])

H3K4me3_length <- dat2$H3K4me3_peak_length
Broad_K4me3_genes<-as.character(dat2$Gene[dat2$H3K4me3_peak_length>4000])
nonBroad_K4me3_genes<-as.character(dat2$Gene[-which(dat2$Gene%in%Broad_K4me3_genes)])


gbdiffmethyl <- read.table("data/genebody_differential_methyl_COSMIC.txt", header=T, sep="\t",fill=TRUE,quote = "")
Cancer_median <- as.numeric(paste(gbdiffmethyl$Cancer_median))
Normal_median <- as.numeric(paste(gbdiffmethyl$Normal_median))
ratio<-Cancer_median/Normal_median

hypermethylated_genebody<-which(ratio > 4)
hypermethylated_genes <- as.character(gbdiffmethyl[hypermethylated_genebody,1])
nonhypermethylated_genes <- as.character(dat2$Gene[-which(anno$Gene%in%hypermethylated_genes)])

quantile_table<-quantile(dat2$LoF_mutations_kb, prob = seq(0, 1, length = 21), type = 5,na.rm=T) #95% percentile
index_high<-which(dat2$LoF_mutations_kb>as.numeric(quantile_table[length(quantile_table)-1]))

high_LoF_genes<-as.character(dat2$Gene[index_high])
nonhigh_LoF_genes<-as.character(dat2$Gene[-which(dat2$Gene%in%high_LoF_genes)])

quantile_table<-quantile(dat2$Missense_mutations_kb, prob = seq(0, 1, length = 21), type = 5,na.rm=T) #95% percentile
index_high<-which(dat2$Missense_mutations_kb>as.numeric(quantile_table[length(quantile_table)-1]))

high_Missense_genes<-as.character(dat2$Gene[index_high])
nonhigh_Missense_genes<-as.character(dat2$Gene[-which(dat2$Gene%in%high_Missense_genes)])

HKG <- read.table("data/HKG_processed.txt", header=T, sep="\t")
HKG_genes<-as.character(HKG[,1])
non_HKG_genes<-as.character(anno$Gene[-which(anno$Gene%in%HKG_genes)])

quantile_table<-quantile(dat2$Missense_o_e_constraint, prob = seq(0, 1, length = 21), type = 5,na.rm=T) #95% percentile
index_high<-which(dat2$Missense_o_e_constraint>as.numeric(quantile_table[length(quantile_table)-1]))

high_z_missense_o_e_genes<-as.character(dat2$Gene[index_high])
nonhigh_z_missense_o_e_genes<-as.character(dat2$Gene[-which(dat2$Gene%in%high_z_missense_o_e_genes)])

quantile_table<-quantile(dat2$LoF_o_e_constraint, prob = seq(0, 1, length = 21), type = 5,na.rm=T) #95% percentile
index_high<-which(dat2$LoF_o_e_constraint>as.numeric(quantile_table[length(quantile_table)-1]))

high_LOF_o_e_constraint_Zscore_genes<-as.character(dat2$Gene[index_high])
nonhigh_LOF_o_e_constraint_Zscore_genes<-as.character(dat2$Gene[-which(dat2$Gene%in%high_LOF_o_e_constraint_Zscore_genes)])

CGC_TSG_genes<-as.character(dat2$Gene[index_TSG])
nonCGC_TSG_genes<-as.character(dat2$Gene[-index_TSG])
CGC_OG_genes<-as.character(dat2$Gene[index_OG])
nonCGC_OG_genes<-as.character(dat2$Gene[-index_OG])
CGC_dual_genes<-as.character(dat2$Gene[index_dual])
nonCGC_dual_genes<-as.character(dat2$Gene[-index_dual])
prioritized_TSG_genes<-as.character(dat2$Gene[index_pTSG])
nonprioritized_TSG_genes<-as.character(dat2$Gene[-index_pTSG])
prioritized_OG_genes<-as.character(dat2$Gene[index_pOG])
nonprioritized_OG_genes<-as.character(dat2$Gene[-index_pOG])
prioritized_dual_genes<-as.character(dat2$Gene[index_pdual])
nonprioritized_dual_genes<-as.character(dat2$Gene[-index_pdual])
NG_genes<-as.character(dat2$Gene[index_NG])
nonNG_genes<-as.character(dat2$Gene[-index_NG])

knownTSG_hub<-as.character(CGC_TSG_genes[which(CGC_TSG_genes%in%hub_genes)])
knownTSG_nonhub<-as.character(CGC_TSG_genes[which(CGC_TSG_genes%in%nonhub_genes)])
nonknownTSG_hub<-as.character(nonCGC_TSG_genes[which(nonCGC_TSG_genes%in%hub_genes)])
nonknownTSG_nonhub<-as.character(nonCGC_TSG_genes[which(nonCGC_TSG_genes%in%nonhub_genes)])
fisher_knownTSG<-matrix(c(length(knownTSG_hub),length(knownTSG_nonhub),length(nonknownTSG_hub), length(nonknownTSG_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("TSG", "nonTSG")))

knownOG_hub<-as.character(CGC_OG_genes[which(CGC_OG_genes%in%hub_genes)])
knownOG_nonhub<-as.character(CGC_OG_genes[which(CGC_OG_genes%in%nonhub_genes)])
nonknownOG_hub<-as.character(nonCGC_OG_genes[which(nonCGC_OG_genes%in%hub_genes)])
nonknownOG_nonhub<-as.character(nonCGC_OG_genes[which(nonCGC_OG_genes%in%nonhub_genes)])
fisher_knownOG<-matrix(c(length(knownOG_hub),length(knownOG_nonhub),length(nonknownOG_hub), length(nonknownOG_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("OG", "nonOG")))

CGCdual_hub<-as.character(CGC_dual_genes[which(CGC_dual_genes%in%hub_genes)])
CGCdual_nonhub<-as.character(CGC_dual_genes[which(CGC_dual_genes%in%nonhub_genes)])
nonCGCdual_hub<-as.character(nonCGC_dual_genes[which(nonCGC_dual_genes%in%hub_genes)])
nonCGCdual_nonhub<-as.character(nonCGC_dual_genes[which(nonCGC_dual_genes%in%nonhub_genes)])
fisher_CGCdual<-matrix(c(length(CGCdual_hub),length(CGCdual_nonhub),length(nonCGCdual_hub), length(nonCGCdual_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("dual", "nondual")))

prioritizedTSG_hub<-as.character(prioritized_TSG_genes[which(prioritized_TSG_genes%in%hub_genes)])
prioritizedTSG_nonhub<-as.character(prioritized_TSG_genes[which(prioritized_TSG_genes%in%nonhub_genes)])
nonprioritizedTSG_hub<-as.character(nonprioritized_TSG_genes[which(nonprioritized_TSG_genes%in%hub_genes)])
nonprioritizedTSG_nonhub<-as.character(nonprioritized_TSG_genes[which(nonprioritized_TSG_genes%in%nonhub_genes)])
fisher_prioritizedTSG<-matrix(c(length(prioritizedTSG_hub),length(prioritizedTSG_nonhub),length(nonprioritizedTSG_hub), length(nonprioritizedTSG_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("prioritized_TSG", "nonprioritized_TSG")))

prioritizedOG_hub<-as.character(prioritized_OG_genes[which(prioritized_OG_genes%in%hub_genes)])
prioritizedOG_nonhub<-as.character(prioritized_OG_genes[which(prioritized_OG_genes%in%nonhub_genes)])
nonprioritizedOG_hub<-as.character(nonprioritized_OG_genes[which(nonprioritized_OG_genes%in%hub_genes)])
nonprioritizedOG_nonhub<-as.character(nonprioritized_OG_genes[which(nonprioritized_OG_genes%in%nonhub_genes)])
fisher_prioritizedOG<-matrix(c(length(prioritizedOG_hub),length(prioritizedOG_nonhub),length(nonprioritizedOG_hub), length(nonprioritizedOG_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("prioritized_OG", "nonprioritized_OG")))

prioritized_dual_hub<-as.character(prioritized_dual_genes[which(prioritized_dual_genes%in%hub_genes)])
prioritized_dual_nonhub<-as.character(prioritized_dual_genes[which(prioritized_dual_genes%in%nonhub_genes)])
nonprioritized_dual_hub<-as.character(nonprioritized_dual_genes[which(nonprioritized_dual_genes%in%hub_genes)])
nonprioritized_dual_nonhub<-as.character(nonprioritized_dual_genes[which(nonprioritized_dual_genes%in%nonhub_genes)])
fisher_prioritized_dual<-matrix(c(length(prioritized_dual_hub),length(prioritized_dual_nonhub),length(nonprioritized_dual_hub), length(nonprioritized_dual_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("prioritized_dual", "nonprioritized_dual")))

NG_hub<-as.character(NG_genes[which(NG_genes%in%hub_genes)])
NG_nonhub<-as.character(NG_genes[which(NG_genes%in%nonhub_genes)])
nonNG_hub<-as.character(nonNG_genes[which(nonNG_genes%in%hub_genes)])
nonNG_nonhub<-as.character(nonNG_genes[which(nonNG_genes%in%nonhub_genes)])
fisher_NG<-matrix(c(length(NG_hub),length(NG_nonhub),length(nonNG_hub), length(nonNG_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("NG", "nonNG")))

essential_hub<-as.character(essential[which(essential%in%hub_genes)])
essential_nonhub<-as.character(essential[which(essential%in%nonhub_genes)])
nonessential_hub<-as.character(non_essential[which(non_essential%in%hub_genes)])
nonessential_nonhub<-as.character(non_essential[which(non_essential%in%nonhub_genes)])
fisher_essential<-matrix(c(length(essential_hub),length(essential_nonhub),length(nonessential_hub), length(nonessential_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("essential", "nonessential")))

canyon_hub<-as.character(canyon_genes[which(canyon_genes%in%hub_genes)])
canyon_nonhub<-as.character(canyon_genes[which(canyon_genes%in%nonhub_genes)])
noncanyon_hub<-as.character(noncanyon_genes[which(noncanyon_genes%in%hub_genes)])
noncanyon_nonhub<-as.character(noncanyon_genes[which(noncanyon_genes%in%nonhub_genes)])
fisher_canyon<-matrix(c(length(canyon_hub),length(canyon_nonhub),length(noncanyon_hub), length(noncanyon_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("canyon", "noncanyon")))

Broad_K4me3_hub<-as.character(Broad_K4me3_genes[which(Broad_K4me3_genes%in%hub_genes)])
Broad_K4me3_nonhub<-as.character(Broad_K4me3_genes[which(Broad_K4me3_genes%in%nonhub_genes)])
nonBroad_K4me3_hub<-as.character(nonBroad_K4me3_genes[which(nonBroad_K4me3_genes%in%hub_genes)])
nonBroad_K4me3_nonhub<-as.character(nonBroad_K4me3_genes[which(nonBroad_K4me3_genes%in%nonhub_genes)])
fisher_Broad_K4me3<-matrix(c(length(Broad_K4me3_hub),length(Broad_K4me3_nonhub),length(nonBroad_K4me3_hub), length(nonBroad_K4me3_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("Broad_K4me3", "nonBroad_K4me3")))

ER_genes_hub<-as.character(ER_genes[which(ER_genes%in%hub_genes)])
ER_genes_nonhub<-as.character(ER_genes[which(ER_genes%in%nonhub_genes)])
nonER_genes_hub<-as.character(non_ER_genes[which(non_ER_genes%in%hub_genes)])
nonER_genes_nonhub<-as.character(non_ER_genes[which(non_ER_genes%in%nonhub_genes)])
fisher_ER_genes<-matrix(c(length(ER_genes_hub),length(ER_genes_nonhub),length(nonER_genes_hub), length(nonER_genes_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("ER_genes", "nonER_genes")))

hypermethylated_genes_hub<-as.character(hypermethylated_genes[which(hypermethylated_genes%in%hub_genes)])
hypermethylated_genes_nonhub<-as.character(hypermethylated_genes[which(hypermethylated_genes%in%nonhub_genes)])
nonhypermethylated_genes_hub<-as.character(nonhypermethylated_genes[which(nonhypermethylated_genes%in%hub_genes)])
nonhypermethylated_genes_nonhub<-as.character(nonhypermethylated_genes[which(nonhypermethylated_genes%in%nonhub_genes)])
fisher_hypermethylated_genes<-matrix(c(length(hypermethylated_genes_hub),length(hypermethylated_genes_nonhub),length(nonhypermethylated_genes_hub), length(nonhypermethylated_genes_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("hypermethylated_genes", "nonhypermethylated_genes")))

high_LoF_genes_hub<-as.character(high_LoF_genes[which(high_LoF_genes%in%hub_genes)])
high_LoF_genes_nonhub<-as.character(high_LoF_genes[which(high_LoF_genes%in%nonhub_genes)])
nonhigh_LoF_genes_hub<-as.character(nonhigh_LoF_genes[which(nonhigh_LoF_genes%in%hub_genes)])
nonhigh_LoF_genes_nonhub<-as.character(nonhigh_LoF_genes[which(nonhigh_LoF_genes%in%nonhub_genes)])
fisher_high_LoF_genes<-matrix(c(length(high_LoF_genes_hub),length(high_LoF_genes_nonhub),length(nonhigh_LoF_genes_hub), length(nonhigh_LoF_genes_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("high_LoF_genes", "nonhigh_LoF_genes")))

high_Missense_genes_hub<-as.character(high_Missense_genes[which(high_Missense_genes%in%hub_genes)])
high_Missense_genes_nonhub<-as.character(high_Missense_genes[which(high_Missense_genes%in%nonhub_genes)])
nonhigh_Missense_genes_hub<-as.character(nonhigh_Missense_genes[which(nonhigh_Missense_genes%in%hub_genes)])
nonhigh_Missense_genes_nonhub<-as.character(nonhigh_Missense_genes[which(nonhigh_Missense_genes%in%nonhub_genes)])
fisher_high_Missense_genes<-matrix(c(length(high_Missense_genes_hub),length(high_Missense_genes_nonhub),length(nonhigh_Missense_genes_hub), length(nonhigh_Missense_genes_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("high_Missense_genes", "nonhigh_Missense_genes")))

HKG_genes_hub<-as.character(HKG_genes[which(HKG_genes%in%hub_genes)])
HKG_genes_nonhub<-as.character(HKG_genes[which(HKG_genes%in%nonhub_genes)])
nonHKG_genes_hub<-as.character(non_HKG_genes[which(non_HKG_genes%in%hub_genes)])
nonHKG_genes_nonhub<-as.character(non_HKG_genes[which(non_HKG_genes%in%nonhub_genes)])
fisher_HKG_genes<-matrix(c(length(HKG_genes_hub),length(HKG_genes_nonhub),length(nonHKG_genes_hub), length(nonHKG_genes_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("HKG_genes", "nonHKG_genes")))

high_z_missense_o_e_genes_hub<-as.character(high_z_missense_o_e_genes[which(high_z_missense_o_e_genes%in%hub_genes)])
high_z_missense_o_e_genes_nonhub<-as.character(high_z_missense_o_e_genes[which(high_z_missense_o_e_genes%in%nonhub_genes)])
nonhigh_z_missense_o_e_genes_hub<-as.character(nonhigh_z_missense_o_e_genes[which(nonhigh_z_missense_o_e_genes%in%hub_genes)])
nonhigh_z_missense_o_e_genes_nonhub<-as.character(nonhigh_z_missense_o_e_genes[which(nonhigh_z_missense_o_e_genes%in%nonhub_genes)])
fisher_high_z_missense_o_e_genes<-matrix(c(length(high_z_missense_o_e_genes_hub),length(high_z_missense_o_e_genes_nonhub),length(nonhigh_z_missense_o_e_genes_hub), length(nonhigh_z_missense_o_e_genes_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("high_z_missense_o_e_genes", "nonhigh_z_missense_o_e_genes")))

high_LOF_o_e_constraint_Zscore_genes_hub<-as.character(high_LOF_o_e_constraint_Zscore_genes[which(high_LOF_o_e_constraint_Zscore_genes%in%hub_genes)])
high_LOF_o_e_constraint_Zscore_genes_nonhub<-as.character(high_LOF_o_e_constraint_Zscore_genes[which(high_LOF_o_e_constraint_Zscore_genes%in%nonhub_genes)])
nonhigh_LOF_o_e_constraint_Zscore_genes_hub<-as.character(nonhigh_LOF_o_e_constraint_Zscore_genes[which(nonhigh_LOF_o_e_constraint_Zscore_genes%in%hub_genes)])
nonhigh_LOF_o_e_constraint_Zscore_genes_nonhub<-as.character(nonhigh_LOF_o_e_constraint_Zscore_genes[which(nonhigh_LOF_o_e_constraint_Zscore_genes%in%nonhub_genes)])
fisher_high_LOF_o_e_constraint_Zscore_genes<-matrix(c(length(high_LOF_o_e_constraint_Zscore_genes_hub),length(high_LOF_o_e_constraint_Zscore_genes_nonhub),length(nonhigh_LOF_o_e_constraint_Zscore_genes_hub), length(nonhigh_LOF_o_e_constraint_Zscore_genes_nonhub)),nrow = 2,dimnames =list(c("hub", "nonhub"),c("high_LOF_o_e_constraint_Zscore_genes", "nonhigh_LOF_o_e_constraint_Zscore_genes")))

pval_knownTSG<-fisher.test(round(fisher_knownTSG*200/length(index_TSG)), alternative = "g")$p.value
pval_knownOG<-fisher.test(round(fisher_knownOG*200/length(index_OG)), alternative = "g")$p.value
pval_CGCdual<-fisher.test(fisher_CGCdual*200/length(index_dual), alternative = "g")$p.value
pval_prioritizedTSG<-fisher.test(round(fisher_prioritizedTSG*200/length(index_pTSG)), alternative = "g")$p.value
pval_prioritizedOG<-fisher.test(round(fisher_prioritizedOG*200/length(index_pOG)), alternative = "g")$p.value
pval_prioritized_dual<-fisher.test(round(fisher_prioritized_dual*200/length(index_pdual)), alternative = "g")$p.value
pval_NG<-fisher.test(round(fisher_NG*200/length(index_NG)), alternative = "g")$p.value
pval_essential<-fisher.test(round(fisher_essential*200/length(essential)), alternative = "g")$p.value
pval_canyon<-fisher.test(round(fisher_canyon*200/length(canyon_genes)), alternative = "g")$p.value
pval_Broad_K4me3<-fisher.test(round(fisher_Broad_K4me3*200/length(Broad_K4me3_genes)), alternative = "g")$p.value
pval_ER_genes<-fisher.test(round(fisher_ER_genes*200/length(ER_genes)), alternative = "g")$p.value
pval_hypermethylated_genes<-fisher.test(round(fisher_hypermethylated_genes*200/length(hypermethylated_genes)), alternative = "g")$p.value
pval_high_LoF_genes<-fisher.test(round(fisher_high_LoF_genes*200/length(high_LoF_genes)), alternative = "g")$p.value
pval_high_Missense_genes<-fisher.test(round(fisher_high_Missense_genes*200/length(high_Missense_genes)), alternative = "g")$p.value
pval_HKG_genes<-fisher.test(round(fisher_HKG_genes*200/length(HKG_genes)), alternative = "g")$p.value
pval_high_z_missense_o_e_genes<-fisher.test(round(fisher_high_z_missense_o_e_genes*200/length(high_z_missense_o_e_genes)), alternative = "g")$p.value
pval_high_LOF_o_e_genes<-fisher.test(round(fisher_high_LOF_o_e_constraint_Zscore_genes*200/length(high_LOF_o_e_constraint_Zscore_genes)), alternative = "g")$p.value

ratio<-length(hub_genes)/length(nonhub_genes)
Gene_number<-c(length(knownTSG_hub)/length(index_TSG),length(knownTSG_nonhub)*ratio/length(index_TSG),length(knownOG_hub)/length(index_OG),length(knownOG_nonhub)*ratio/length(index_OG),length(CGCdual_hub)/length(index_dual),length(CGCdual_nonhub)*ratio/length(index_dual),length(prioritizedTSG_hub)/length(index_pTSG),length(prioritizedTSG_nonhub)*ratio/length(index_pTSG),length(prioritizedOG_hub)/length(index_pOG),length(prioritizedOG_nonhub)*ratio/length(index_pOG),length(prioritized_dual_hub)/length(index_pdual),length(prioritized_dual_nonhub)*ratio/length(index_pdual),length(NG_hub)/length(index_NG),length(NG_nonhub)*ratio/length(index_NG),length(essential_hub)/length(essential),length(essential_nonhub)*ratio/length(essential),length(canyon_hub)/length(canyon_genes),length(canyon_nonhub)*ratio/length(canyon_genes),length(Broad_K4me3_hub)/length(Broad_K4me3_genes),length(Broad_K4me3_nonhub)*ratio/length(Broad_K4me3_genes),
length(ER_genes_hub)/length(ER_genes),length(ER_genes_nonhub)*ratio/length(ER_genes),length(hypermethylated_genes_hub)/length(hypermethylated_genes),length(hypermethylated_genes_nonhub)*ratio/length(hypermethylated_genes),length(high_LoF_genes_hub)/length(high_LoF_genes),length(high_LoF_genes_nonhub)*ratio/length(high_LoF_genes),length(high_Missense_genes_hub)/length(high_Missense_genes),length(high_Missense_genes_nonhub)*ratio/length(high_Missense_genes),length(HKG_genes_hub)/length(HKG_genes),length(HKG_genes_nonhub)*ratio/length(HKG_genes),length(high_z_missense_o_e_genes_hub)/length(high_z_missense_o_e_genes),length(high_z_missense_o_e_genes_nonhub)*ratio/length(high_z_missense_o_e_genes),length(high_LOF_o_e_constraint_Zscore_genes_hub)/length(high_LOF_o_e_constraint_Zscore_genes),length(high_LOF_o_e_constraint_Zscore_genes_nonhub)*ratio/length(high_LOF_o_e_constraint_Zscore_genes))*1000

#index_large<-which(Gene_number>150)
#Gene_number <- as.numeric(Gene_number)
#Gene_number[index_large]<-(Gene_number[index_large]-150)/2+150

df1<-read.table(header=T,text='stage	genetype
Hub	CGC-TSG
non-Hub	CGC-TSG
Hub	CGC-OG
non-Hub	CGC-OG
Hub	CGC-dual
non-Hub	CGC-dual
Hub	"Novel DORGE-TSG"
non-Hub	"Novel DORGE-TSG"
Hub	"Novel DORGE-OG"
non-Hub	"Novel DORGE-OG"
Hub	"Novel DORGE-dual"
non-Hub	"Novel DORGE-dual"
Hub	NG
non-Hub	NG
');
df1$genetype<-factor(df1$genetype,levels=c("NG", "CGC-dual", "CGC-OG", "CGC-TSG", "Novel DORGE-dual", "Novel DORGE-OG", "Novel DORGE-TSG"))
df_2<-cbind(df1,Gene_number=Gene_number[1:nrow(df1)])


df<-read.table(header=T,text='stage	genetype
Hub	Essential
non-Hub	Essential
Hub	"Gene-body canyon"
non-Hub	"Gene-body canyon"
Hub	"Broad H3K4me3"
non-Hub	"Broad H3K4me3"
Hub	ER
non-Hub	ER
Hub	"Gene-body\n hypermethylation"
non-Hub	"Gene-body\n hypermethylation"
Hub	"High LoF/kb"
non-Hub	"High LoF/kb"
Hub	"High missense/kb"
non-Hub	"High missense/kb"
Hub	HKG
non-Hub	HKG
Hub	"High missense\n o/e constraint"
non-Hub	"High missense\n o/e constraint"
Hub	"High LoF o/e constraint"
non-Hub	"High LoF o/e constraint"
');
df$genetype<-factor(df$genetype,levels=c("Essential", "HKG", "Broad H3K4me3"
, "ER", "Gene-body\n hypermethylation", "Gene-body canyon", "High LoF/kb","High LoF o/e constraint", "High missense/kb"
, "High missense\n o/e constraint"))
df_3<-cbind(df,Gene_number=Gene_number[(nrow(df1)+1):length(Gene_number)])

pdf("Raw_figures/Figure_5C_BioGRID_network_hub_enrichment.pdf", family="ArialMT", width=4.19814*1.12, height=3.0)
p<-ggplot(data=df_2, aes(x=genetype, y=Gene_number, fill=stage))+geom_bar(stat="identity", position=position_dodge())+theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.position="top")+ theme_bw() + theme(axis.text.y = element_text(hjust = 1, colour = "black"),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),plot.title = element_text(hjust = 0.5),legend.title = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual("legend", values = c("Hub" = "#E69F00", "non-Hub" = "grey"))+ labs(x = "", y = "Gene number")+ggtitle("PPI network")+ylim(0,350)
p+geom_signif(y_position=c(max(df_2[13,3],df_2[14,3])+10, max(df_2[5,3],df_2[6,3])+10), xmin=c(0.8,1.8), xmax=c(1.2,2.2),annotation=c(formatC(pval_NG, digits = 2), formatC(pval_CGCdual, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)+geom_signif(y_position=c(max(df_2[3,3],df_2[4,3])+15, max(df_2[1,3],df_2[2,3])+10), xmin=c(2.8,3.8), xmax=c(3.2,4.2),annotation=c(formatC(pval_knownOG, format = "e", digits = 2), formatC(pval_knownTSG, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)+geom_signif(y_position=c(max(df_2[11,3],df_2[12,3])+10, max(df_2[9,3],df_2[10,3])+10), xmin=c(4.8,5.8), xmax=c(5.2,6.2),annotation=c(formatC(pval_prioritized_dual, format = "e", digits = 2), formatC(pval_prioritizedOG, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)+geom_signif(y_position=c(max(df_2[7,3],df_2[8,3])+10), xmin=c(6.8), xmax=c(7.2),annotation=c(formatC(pval_prioritizedTSG, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)
garbage <- dev.off()

pdf("Raw_figures/Figure_5D_BioGRID_network_hub_enrichment.pdf", family="ArialMT", width=5.2196*1.12, height=3)
p<-ggplot(data=df_3, aes(x=genetype, y=Gene_number, fill=stage))+geom_bar(stat="identity", position=position_dodge())+theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.position="top")+ theme_bw() + theme(axis.text.y = element_text(hjust = 1, colour = "black"),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),plot.title = element_text(hjust = 0.5),legend.title = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual("legend", values = c("Hub" = "#E69F00", "non-Hub" = "grey"))+ labs(x = "", y = "Gene number")+ggtitle("PPI network")+ylim(0,350)
p+geom_signif(y_position=c(max(df_3[1,3],df_3[2,3])+10,max(df_3[15,3],df_3[16,3])+10), xmin=c(0.8,1.8), xmax=c(1.2,2.2),annotation=c(formatC(pval_essential, format = "e", digits = 2), formatC(pval_HKG_genes, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)+geom_signif(y_position=c(max(df_3[5,3],df_3[6,3])+10, max(df_3[7,3],df_3[8,3])+10), xmin=c(2.8,3.8), xmax=c(3.2,4.2),annotation=c(formatC(pval_Broad_K4me3, format = "e", digits = 2), formatC(pval_ER_genes, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)+geom_signif(y_position=c(max(df_3[9,3],df_3[10,3])+10, max(df_3[3,3],df_3[4,3])+10), xmin=c(4.8,5.8), xmax=c(5.2,6.2),annotation=c(formatC(pval_hypermethylated_genes, digits = 2), formatC(pval_canyon, digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)+geom_signif(y_position=c(max(df_3[11,3],df_3[12,3])+10, max(df_3[19,3],df_3[20,3])+10), xmin=c(6.8,7.8), xmax=c(7.2,8.2),annotation=c(formatC(pval_high_LoF_genes, format = "e", digits = 2), formatC(pval_high_LOF_o_e_genes, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)+geom_signif(y_position=c(max(df_3[13,3],df_3[14,3])+10, max(df_3[17,3],df_3[18,3])+10), xmin=c(8.8,9.8), xmax=c(9.2,10.2),annotation=c(formatC(pval_high_Missense_genes, format = "e", digits = 2), formatC(pval_high_z_missense_o_e_genes, format = "e", digits = 2)), tip_length=0.02,size=0.3,textsize =1.7)

garbage <- dev.off()