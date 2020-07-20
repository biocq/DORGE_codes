###################
##### Heatmap #####
###################

###### The gene enrichment analyses for different gene categories to independently evaluate DORGE prediction.

###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/Phyletic/ Eukaryota_4.txt, Metazoa_3.txt, Chordata_2.txt, Mammalia_1.txt (Processed Phyletic age gene lists (gene age) downloaded from Online GEne Essentiality (OGEE) database, only human records were analyzed)
###### Output: Figure_3G_Phyletic_enrichment.pdf: Enrichment of different Phyletic age genes in different gene categories


options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()
pkgs <-  c("ComplexHeatmap")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
}

pkgs <-  c("ggpubr","circlize")


if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("ggpubr"))
suppressMessages(library("circlize"))
suppressMessages(library("ComplexHeatmap"))

TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

##################### Figure 3G: Phyletic age gene set enrichment #####################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_3");

anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","TSG_all","OG_all","NG")]

index_NG<-which(allgene$NG=="1")
index_pTSG<-which(allgene$TSG_probability>TSG_threshold & allgene$OG_probability< OG_threshold & allgene$TSG_all!="1")
index_pOG<-which(allgene$OG_probability> OG_threshold & allgene$TSG_probability<TSG_threshold & allgene$OG_all!="1")
index_TSG<-which(allgene$TSG_all=="1" & allgene$OG_all!="1")
index_OG<-which(allgene$OG_all=="1" & allgene$TSG_all!="1")

known_tsg_genes<-as.character(allgene$Gene[index_TSG])
known_og_genes<-as.character(allgene$Gene[index_OG])
new_predicted_tsgs<-as.character(allgene$Gene[index_pTSG])
new_predicted_ogs<-as.character(allgene$Gene[index_pOG])
NG<-as.character(allgene$Gene[anno[,"NG"]==1])

nonknownTSG<-setdiff(allgene[,1],known_tsg_genes)
nonknownOG<-setdiff(allgene[,1],known_og_genes)
nonPredictedTSG<-setdiff(allgene[,1],new_predicted_tsgs)
nonPredictedOG<-setdiff(allgene[,1],new_predicted_ogs)
nonNG<-setdiff(allgene[,1],NG)

#Eukaryota

Eukaryota <- read.table("data/Phyletic/Eukaryota_4.txt", header=F, sep="\t")
Eukaryota_genes <- as.character(allgene$Gene[which(as.character(allgene$Gene)%in%Eukaryota$V1)])
nonEukaryota <- as.character(allgene$Gene[-which(as.character(allgene$Gene)%in%Eukaryota$V1)])

#Predicted TSG
predictedTSG_Eukaryota<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%Eukaryota_genes)])
predictedTSG_nonEukaryota<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%nonEukaryota)])
nonPredictedTSG_Eukaryota<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%Eukaryota_genes)])
nonPredictedTSG_nonEukaryota<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%nonEukaryota)])
fisher<-matrix(c(length(predictedTSG_Eukaryota),length(predictedTSG_nonEukaryota),length(nonPredictedTSG_Eukaryota), length(nonPredictedTSG_nonEukaryota))*200/length(new_predicted_tsgs),nrow = 2,dimnames =list(c("Eukaryota", "nonEukaryota"),c("TSG", "nonTSG")))
pval_predictedTSG_Eukaryota<-fisher.test(round(fisher), alternative = "g")$p.value

#known TSG

knownTSG_Eukaryota<-as.character(known_tsg_genes[which(known_tsg_genes%in%Eukaryota_genes)])
knownTSG_nonEukaryota<-as.character(known_tsg_genes[which(known_tsg_genes%in%nonEukaryota)])
nonknownTSG_Eukaryota<-as.character(nonknownTSG[which(nonknownTSG%in%Eukaryota_genes)])
nonknownTSG_nonEukaryota<-as.character(nonknownTSG[which(nonknownTSG%in%nonEukaryota)])
fisher<-matrix(c(length(knownTSG_Eukaryota),length(knownTSG_nonEukaryota),length(nonknownTSG_Eukaryota), length(nonknownTSG_nonEukaryota))*200/length(new_predicted_ogs),nrow = 2,dimnames =list(c("Eukaryota", "nonEukaryota"),c("TSG", "nonTSG")))
pval_knownTSG_Eukaryota<-fisher.test(round(fisher), alternative = "g")$p.value

#Novel OG

predictedOG_Eukaryota<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%Eukaryota_genes)])
predictedOG_nonEukaryota<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%nonEukaryota)])
nonPredictedOG_Eukaryota<-nonPredictedOG[which(nonPredictedOG%in%Eukaryota_genes)]
nonPredictedOG_nonEukaryota<-as.character(nonPredictedOG[which(nonPredictedOG%in%nonEukaryota)])
fisher<-matrix(c(length(predictedOG_Eukaryota),length(predictedOG_nonEukaryota),length(nonPredictedOG_Eukaryota), length(nonPredictedOG_nonEukaryota))*200/length(new_predicted_ogs),nrow = 2,dimnames =list(c("Eukaryota", "nonEukaryota"),c("OG", "nonOG")))
pval_predictedOG_Eukaryota<-fisher.test(round(fisher), alternative = "g")$p.value

#Known OG

knownOG_Eukaryota<-as.character(known_og_genes[which(known_og_genes%in%Eukaryota_genes)])
knownOG_nonEukaryota<-as.character(known_og_genes[which(known_og_genes%in%nonEukaryota)])
nonknownOG_Eukaryota<-as.character(nonknownOG[which(nonknownOG%in%Eukaryota_genes)])
nonknownOG_nonEukaryota<-as.character(nonknownOG[which(nonknownOG%in%nonEukaryota)])
fisher<-matrix(c(length(knownOG_Eukaryota),length(knownOG_nonEukaryota),length(nonknownOG_Eukaryota), length(nonknownOG_nonEukaryota))*200/length(known_og_genes),nrow = 2,dimnames =list(c("Eukaryota", "nonEukaryota"),c("OG", "nonOG")))
pval_knownOG_Eukaryota<-fisher.test(round(fisher), alternative = "g")$p.value

#NG

NG_Eukaryota<-as.character(NG[which(NG%in%Eukaryota_genes)])
NG_nonEukaryota<-as.character(NG[which(NG%in%nonEukaryota)])
nonNG_Eukaryota<-as.character(nonNG[which(nonNG%in%Eukaryota_genes)])
nonNG_nonEukaryota<-as.character(nonNG[which(nonNG%in%nonEukaryota)])
fisher<-matrix(c(length(NG_Eukaryota),length(NG_nonEukaryota),length(nonNG_Eukaryota), length(NG))*200/length(NG),nrow = 2,dimnames =list(c("Eukaryota", "nonEukaryota"),c("NG", "nonNG")))
pval_NG_Eukaryota<-fisher.test(round(fisher), alternative = "g")$p.value


#Metazoa
Metazoa <- read.table("data/Phyletic/Metazoa_3.txt", header=F, sep="\t")
Metazoa_genes <- allgene$Gene[which(as.character(allgene$Gene)%in%Metazoa$V1)]
nonMetazoa <- allgene$Gene[-which(as.character(allgene$Gene)%in%Metazoa$V1)]

#Predicted TSG
predictedTSG_Metazoa<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%Metazoa_genes)])
predictedTSG_nonMetazoa<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%nonMetazoa)])
nonPredictedTSG_Metazoa<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%Metazoa_genes)])
nonPredictedTSG_nonMetazoa<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%nonMetazoa)])
fisher<-matrix(c(length(predictedTSG_Metazoa),length(predictedTSG_nonMetazoa),length(nonPredictedTSG_Metazoa), length(nonPredictedTSG_nonMetazoa))*200/length(new_predicted_tsgs),nrow = 2,dimnames =list(c("Metazoa", "nonMetazoa"),c("TSG", "nonTSG")))
pval_predictedTSG_Metazoa<-fisher.test(round(fisher), alternative = "g")$p.value

#known TSG

knownTSG_Metazoa<-as.character(known_tsg_genes[which(known_tsg_genes%in%Metazoa_genes)])
knownTSG_nonMetazoa<-as.character(known_tsg_genes[which(known_tsg_genes%in%nonMetazoa)])
nonknownTSG_Metazoa<-as.character(nonknownTSG[which(nonknownTSG%in%Metazoa_genes)])
nonknownTSG_nonMetazoa<-as.character(nonknownTSG[which(nonknownTSG%in%nonMetazoa)])
fisher<-matrix(c(length(knownTSG_Metazoa),length(knownTSG_nonMetazoa),length(nonknownTSG_Metazoa), length(nonknownTSG_nonMetazoa))*200/length(known_tsg_genes),nrow = 2,dimnames =list(c("Metazoa", "nonMetazoa"),c("TSG", "nonTSG")))
pval_knownTSG_Metazoa<-fisher.test(round(fisher), alternative = "g")$p.value

#Novel OG

predictedOG_Metazoa<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%Metazoa_genes)])
predictedOG_nonMetazoa<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%nonMetazoa)])
nonPredictedOG_Metazoa<-as.character(nonPredictedOG[which(nonPredictedOG%in%Metazoa_genes)])
nonPredictedOG_nonMetazoa<-as.character(nonPredictedOG[which(nonPredictedOG%in%nonMetazoa)])
fisher<-matrix(c(length(predictedOG_Metazoa),length(predictedOG_nonMetazoa),length(nonPredictedOG_Metazoa), length(nonPredictedOG_nonMetazoa))*200/length(new_predicted_ogs),nrow = 2,dimnames =list(c("Metazoa", "nonMetazoa"),c("OG", "nonOG")))
pval_predictedOG_Metazoa<-fisher.test(round(fisher), alternative = "g")$p.value

#Known OG

knownOG_Metazoa<-as.character(known_og_genes[which(known_og_genes%in%Metazoa_genes)])
knownOG_nonMetazoa<-as.character(known_og_genes[which(known_og_genes%in%nonMetazoa)])
nonknownOG_Metazoa<-as.character(nonknownOG[which(nonknownOG%in%Metazoa_genes)])
nonknownOG_nonMetazoa<-as.character(nonknownOG[which(nonknownOG%in%nonMetazoa)])
fisher<-matrix(c(length(knownOG_Metazoa),length(knownOG_nonMetazoa),length(nonknownOG_Metazoa), length(nonknownOG_nonMetazoa))*200/length(known_og_genes),nrow = 2,dimnames =list(c("Metazoa", "nonMetazoa"),c("OG", "nonOG")))
pval_knownOG_Metazoa<-fisher.test(round(fisher), alternative = "g")$p.value

#NG

NG_Metazoa<-as.character(NG[which(NG%in%Metazoa_genes)])
NG_nonMetazoa<-as.character(NG[which(NG%in%nonMetazoa)])
nonNG_Metazoa<-as.character(nonNG[which(nonNG%in%Metazoa_genes)])
nonNG_nonMetazoa<-as.character(nonNG[which(nonNG%in%nonMetazoa)])
fisher<-matrix(c(length(NG_Metazoa),length(NG_nonMetazoa),length(nonNG_Metazoa), length(nonNG_nonMetazoa))*200/length(NG),nrow = 2,dimnames =list(c("Metazoa", "nonMetazoa"),c("NG", "nonNG")))
pval_NG_Metazoa<-fisher.test(round(fisher), alternative = "g")$p.value

#Chordata

Chordata <- read.table("data/Phyletic/Chordata_2.txt", header=F, sep="\t")

#novel TSG
Chordata_genes <- allgene$Gene[which(as.character(allgene$Gene)%in%Chordata$V1)]
nonChordata <- allgene$Gene[-which(as.character(allgene$Gene)%in%Chordata$V1)]

#Predicted TSG
predictedTSG_Chordata<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%Chordata_genes)])
predictedTSG_nonChordata<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%nonChordata)])
nonPredictedTSG_Chordata<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%Chordata_genes)])
nonPredictedTSG_nonChordata<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%nonChordata)])
fisher<-matrix(c(length(predictedTSG_Chordata),length(predictedTSG_nonChordata),length(nonPredictedTSG_Chordata), length(nonPredictedTSG_nonChordata))*200/length(new_predicted_tsgs),nrow = 2,dimnames =list(c("Chordata", "nonChordata"),c("TSG", "nonTSG")))
pval_predictedTSG_Chordata<-fisher.test(round(fisher), alternative = "g")$p.value

#known TSG

knownTSG_Chordata<-as.character(known_tsg_genes[which(known_tsg_genes%in%Chordata_genes)])
knownTSG_nonChordata<-as.character(known_tsg_genes[which(known_tsg_genes%in%nonChordata)])
nonknownTSG_Chordata<-as.character(nonknownTSG[which(nonknownTSG%in%Chordata_genes)])
nonknownTSG_nonChordata<-as.character(nonknownTSG[which(nonknownTSG%in%nonChordata)])
fisher<-matrix(c(length(knownTSG_Chordata),length(knownTSG_nonChordata),length(nonknownTSG_Chordata), length(nonknownTSG_nonChordata))*200/length(known_tsg_genes),nrow = 2,dimnames =list(c("Chordata", "nonChordata"),c("TSG", "nonTSG")))
pval_knownTSG_Chordata<-fisher.test(round(fisher), alternative = "g")$p.value

#Novel OG

predictedOG_Chordata<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%Chordata_genes)])
predictedOG_nonChordata<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%nonChordata)])
nonPredictedOG_Chordata<-as.character(nonPredictedOG[which(nonPredictedOG%in%Chordata_genes)])
nonPredictedOG_nonChordata<-as.character(nonPredictedOG[which(nonPredictedOG%in%nonChordata)])
fisher<-matrix(c(length(predictedOG_Chordata),length(predictedOG_nonChordata),length(nonPredictedOG_Chordata), length(nonPredictedOG_nonChordata))*200/length(new_predicted_ogs),nrow = 2,dimnames =list(c("Chordata", "nonChordata"),c("OG", "nonOG")))
pval_predictedOG_Chordata<-fisher.test(round(fisher), alternative = "g")$p.value

#Known OG

knownOG_Chordata<-as.character(known_og_genes[which(known_og_genes%in%Chordata_genes)])
knownOG_nonChordata<-as.character(known_og_genes[which(known_og_genes%in%nonChordata)])
nonknownOG_Chordata<-as.character(nonknownOG[which(nonknownOG%in%Chordata_genes)])
nonknownOG_nonChordata<-as.character(nonknownOG[which(nonknownOG%in%nonChordata)])
fisher<-matrix(c(length(knownOG_Chordata),length(knownOG_nonChordata),length(nonknownOG_Chordata), length(nonknownOG_nonChordata))*200/length(known_og_genes),nrow = 2,dimnames =list(c("Chordata", "nonChordata"),c("OG", "nonOG")))
pval_knownOG_Chordata<-fisher.test(round(fisher), alternative = "g")$p.value

#NG

NG_Chordata<-as.character(NG[which(NG%in%Chordata_genes)])
NG_nonChordata<-as.character(NG[which(NG%in%nonChordata)])
nonNG_Chordata<-as.character(nonNG[which(nonNG%in%Chordata_genes)])
nonNG_nonChordata<-as.character(nonNG[which(nonNG%in%nonChordata)])
fisher<-matrix(c(length(NG_Chordata),length(NG_nonChordata),length(nonNG_Chordata), length(nonNG_nonChordata))*200/length(NG),nrow = 2,dimnames =list(c("Chordata", "nonChordata"),c("NG", "nonNG")))
pval_NG_Chordata<-fisher.test(round(fisher), alternative = "g")$p.value

#Mammalia

Mammalia <- read.table("data/Phyletic/Mammalia_1.txt", header=F, sep="\t")
Mammalia_genes <- allgene$Gene[which(as.character(allgene$Gene)%in%Mammalia$V1)]
nonMammalia <- allgene$Gene[-which(as.character(allgene$Gene)%in%Mammalia$V1)]

#Predicted TSG
predictedTSG_Mammalia<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%Mammalia_genes)])
predictedTSG_nonMammalia<-as.character(new_predicted_tsgs[which(new_predicted_tsgs%in%nonMammalia)])
nonPredictedTSG_Mammalia<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%Mammalia_genes)])
nonPredictedTSG_nonMammalia<-as.character(nonPredictedTSG[which(nonPredictedTSG%in%nonMammalia)])
fisher<-matrix(c(length(predictedTSG_Mammalia),length(predictedTSG_nonMammalia),length(nonPredictedTSG_Mammalia), length(nonPredictedTSG_nonMammalia))*200/length(new_predicted_tsgs),nrow = 2,dimnames =list(c("Mammalia", "nonMammalia"),c("TSG", "nonTSG")))
pval_predictedTSG_Mammalia<-fisher.test(round(fisher), alternative = "g")$p.value

#known TSG

knownTSG_Mammalia<-as.character(known_tsg_genes[which(known_tsg_genes%in%Mammalia_genes)])
knownTSG_nonMammalia<-as.character(known_tsg_genes[which(known_tsg_genes%in%nonMammalia)])
nonknownTSG_Mammalia<-as.character(nonknownTSG[which(nonknownTSG%in%Mammalia_genes)])
nonknownTSG_nonMammalia<-as.character(nonknownTSG[which(nonknownTSG%in%nonMammalia)])
fisher<-matrix(c(length(knownTSG_Mammalia),length(knownTSG_nonMammalia),length(nonknownTSG_Mammalia), length(nonknownTSG_nonMammalia))*200/length(known_tsg_genes),nrow = 2,dimnames =list(c("Mammalia", "nonMammalia"),c("TSG", "nonTSG")))
pval_knownTSG_Mammalia<-fisher.test(round(fisher), alternative = "g")$p.value

#Novel OG

predictedOG_Mammalia<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%Mammalia_genes)])
predictedOG_nonMammalia<-as.character(new_predicted_ogs[which(new_predicted_ogs%in%nonMammalia)])
nonPredictedOG_Mammalia<-as.character(nonPredictedOG[which(nonPredictedOG%in%Mammalia_genes)])
nonPredictedOG_nonMammalia<-as.character(nonPredictedOG[which(nonPredictedOG%in%nonMammalia)])
fisher<-matrix(c(length(predictedOG_Mammalia),length(predictedOG_nonMammalia),length(nonPredictedOG_Mammalia), length(nonPredictedOG_nonMammalia))*200/length(new_predicted_ogs),nrow = 2,dimnames =list(c("Mammalia", "nonMammalia"),c("OG", "nonOG")))
pval_predictedOG_Mammalia<-fisher.test(round(fisher), alternative = "g")$p.value

#Known OG

knownOG_Mammalia<-as.character(known_og_genes[which(known_og_genes%in%Mammalia_genes)])
knownOG_nonMammalia<-as.character(known_og_genes[which(known_og_genes%in%nonMammalia)])
nonknownOG_Mammalia<-as.character(nonknownOG[which(nonknownOG%in%Mammalia_genes)])
nonknownOG_nonMammalia<-as.character(nonknownOG[which(nonknownOG%in%nonMammalia)])
fisher<-matrix(c(length(knownOG_Mammalia),length(knownOG_nonMammalia),length(nonknownOG_Mammalia), length(nonknownOG_nonMammalia))*200/length(known_og_genes),nrow = 2,dimnames =list(c("Mammalia", "nonMammalia"),c("OG", "nonOG")))
pval_knownOG_Mammalia<-fisher.test(round(fisher), alternative = "g")$p.value

#NG

NG_Mammalia<-as.character(NG[which(NG%in%Mammalia_genes)])
NG_nonMammalia<-as.character(NG[which(NG%in%nonMammalia)])
nonNG_Mammalia<-as.character(nonNG[which(nonNG%in%Mammalia_genes)])
nonNG_nonMammalia<-as.character(nonNG[which(nonNG%in%nonMammalia)])
fisher<-matrix(c(length(NG_Mammalia),length(NG_nonMammalia),length(nonNG_Mammalia), length(nonNG_nonMammalia))*200/length(NG),nrow = 2,dimnames =list(c("Mammalia", "nonMammalia"),c("NG", "nonNG")))
pval_NG_Mammalia<-fisher.test(round(fisher), alternative = "g")$p.value

############### heatmap #################

TSG_pvals<-matrix(c(pval_predictedTSG_Eukaryota,pval_knownTSG_Eukaryota,pval_NG_Eukaryota,pval_predictedTSG_Metazoa,pval_knownTSG_Metazoa,pval_NG_Metazoa,pval_predictedTSG_Chordata,pval_knownTSG_Chordata,pval_NG_Chordata,pval_predictedTSG_Mammalia,pval_knownTSG_Mammalia,pval_NG_Mammalia),nrow=3,ncol=4,byrow = F)
rownames(TSG_pvals)<-c("Novel DORGE-TSG","CGC-TSG","NG")
colnames(TSG_pvals)<-c("Eukaryota","Metazoa","Chordata","Mammalia")
OG_pvals<-matrix(c(pval_predictedOG_Eukaryota,pval_knownOG_Eukaryota,pval_NG_Eukaryota,pval_predictedOG_Metazoa,pval_knownOG_Metazoa,pval_NG_Metazoa,pval_predictedOG_Chordata,pval_knownOG_Chordata,pval_NG_Chordata,pval_predictedOG_Mammalia,pval_knownOG_Mammalia,pval_NG_Mammalia),nrow=3,ncol=4,byrow = F)
rownames(OG_pvals)<-c("Novel DORGE-OG","CGC-OG","NG")
colnames(OG_pvals)<-c("Eukaryota","Metazoa","Chordata","Mammalia")
allpvals<-rbind(TSG_pvals,OG_pvals)
allpvals2 <- t(allpvals)
allpvals2 <- as.data.frame(allpvals2)
log10_pvals<- -log10(allpvals2[,1:5])
log10_pvals[log10_pvals>4]<-4
cn = colnames(log10_pvals)
pdf("Raw_figures/Figure_3G_Phyletic_enrichment.pdf", family="ArialMT", width=6, height=3)
h1<-Heatmap(as.matrix(log10_pvals), row_order =c("Eukaryota","Metazoa","Chordata","Mammalia"), column_order =c("CGC-OG","CGC-TSG","NG","Novel DORGE-OG","Novel DORGE-TSG"),show_column_names=F, col = colorRamp2(c(4, 2,1,0.1), c("red", "darksalmon","deepskyblue","blue")), border = TRUE,name="-log10 P-value",show_column_dend = F,show_row_dend = F,heatmap_height = unit(7,"cm"),heatmap_width = unit(5,"cm"),heatmap_legend_param = list(title_position = "topcenter",color_bar = "continuous",legend_direction = "horizontal",legend_width = unit(3, "cm")),show_row_names = T,bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 35, location = unit(1, "npc"), just = "right"),annotation_height = max_text_width(cn)))
h1
garbage <- dev.off()