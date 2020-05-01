########################
##### Scatterplots #####
########################

###### Scatter plot of H3K4me3 peak length and gene-body differential methylation for CGC genes and DORGE-predicted novel genes

###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: sampled_NGs.txt: Sampled NGs for reproducing purpose

###### Output: Figure_3F_H3K4me3_differential_methylation_scatter_for_DORGE_novel_genes.pdf: Scatter plot of H3K4me3 peak length and gene-body differential methylation for DORGE-predicted novel genes
###### Output: Figure_S3E_H3K4me3_differential_methylation_scatter_for_CGC_genes.pdf: Scatter plot of H3K4me3 peak length and gene-body differential methylation for CGC genes
options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggpubr","scales","dplyr","circlize")


if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}


suppressMessages(library("circlize"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggpubr"))
suppressMessages(library("scales"))

TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005

######### Figure 3F: Scatter plot of H3K4me3 peak length and gene-body differential methylation for DORGE-predicted novel genes #########
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

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_3");

all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","H3K4me3_peak_length","Gene_body_hypermethylation_in_cancer")]

lowdiff<-which(all_feature$Gene_body_hypermethylation_in_cancer< -0.3)
all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]<-(all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]+0.3)/3-0.3
anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
TSG_CGC<-anno[,"TSG_all"]
OG_CGC<-anno[,"OG_all"]
NG<-anno[,"NG"]
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1" #,"PTEN","FOXP1"
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])
allgene<-cbind(allgene,all_feature[,c("H3K4me3_peak_length","Gene_body_hypermethylation_in_cancer")])

allgene$Gene_body_hypermethylation_in_cancer[is.na(allgene$Gene_body_hypermethylation_in_cancer)] <- 0
allgene$H3K4me3_peak_length[is.na(allgene$H3K4me3_peak_length)] <- 0
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>5000)
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-5000)/3+5000

index<-rep("",nrow(allgene))
index_NG<-which(NG=="1")
index_pTSGOG<-which((allgene$TSG_probability>TSG_threshold & OG_CGC=="1" & TSG_CGC!="1")|(allgene$OG_probability> OG_threshold & TSG_CGC=="1" & OG_CGC!="1")| (allgene$TSG_probability>TSG_threshold & allgene$OG_probability> OG_threshold & OG_CGC!="1" & TSG_CGC!="1"))


index_pOG<-c(which(allgene$OG_probability>OG_threshold & allgene$TSG_probability<TSG_threshold & OG_CGC!="1"))
index_pTSG<-which(allgene$TSG_probability>TSG_threshold & allgene$OG_probability<OG_threshold & TSG_CGC!="1")

sampled_index<-sample(which(NG=="1"),2000) # & allgene$OG_probability<OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")

DORGE_novel_index<-rep("",nrow(allgene))
DORGE_novel_index[sampled_index$NG_sampled]<-"NG"
DORGE_novel_index[index_pTSG]<-"Novel DORGE-TSG"
DORGE_novel_index[index_pOG]<-"Novel DORGE-OG"
DORGE_novel_index[index_pTSGOG]<-"Novel DORGE-dual"
DORGE_novel_index <- factor(DORGE_novel_index,levels=c("Novel DORGE-dual", "Novel DORGE-OG", "Novel DORGE-TSG", "NG",""))

allgene1<-cbind(allgene,DORGE_novel_index)
allgene1$DORGE_novel_index <- as.factor(allgene1$DORGE_novel_index)
allgene2<-allgene1[allgene1$DORGE_novel_index!="",]

pdf("Raw_figures/Figure_3F_H3K4me3_differential_methylation_scatter_for_DORGE_novel_genes.pdf", family="ArialMT", width=5, height=4.6179)#Figure 3D
p<-ggscatter(allgene2, x = "H3K4me3_peak_length", y = "Gene_body_hypermethylation_in_cancer",xlab="H3K4me3 peak length", ylab="Gene-body hypermethylation",show.legend.text=F,color = "DORGE_novel_index",size = 1.6,font.label = c(9,"italic"),rug=F,label = "Core",shape=16,label.select = list(criteria = "`x` > 4000 | `y`>0.3"),repel=T,alpha=0.9,palette=c("darkorange2","#4DAF4A","#984EA3","#CCCCCC"),alpha=0.7)+  scale_y_continuous(trans = modulus_trans(8),breaks = c(-0.4,-0.3,0,0.3,0.8),labels = c("-0.6", "-0.3", "0", "0.3", "0.8"))+ scale_x_continuous(trans = modulus_trans(1.6),breaks = c(0,2500,4000,5000,6000),labels = c("0", "2500", "4000", "5000", "8000"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
ggpar(p,legend.title = "",legend=c(0.77,0.82),font.legend = c(7, "plain", "black"))+ rremove("legend.title")
garbage <- dev.off()

######### Figure S3E: Scatter plot of H3K4me3 peak length and gene-body differential methylation for CGC genes #########
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1" #,"PTEN","FOXP1"
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])
allgene<-cbind(allgene,all_feature[,c("H3K4me3_peak_length","Gene_body_hypermethylation_in_cancer")])

allgene$Gene_body_hypermethylation_in_cancer[is.na(allgene$Gene_body_hypermethylation_in_cancer)] <- 0
allgene$H3K4me3_peak_length[is.na(allgene$H3K4me3_peak_length)] <- 0
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>5000)
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-5000)/3+5000

index<-rep("",nrow(allgene))
index_TSG<-which((TSG_CGC=="1" & OG_CGC!="1"))
index_OG<-which((OG_CGC=="1" & TSG_CGC!="1"))
index_TSGOG<-which((OG_CGC=="1" & TSG_CGC=="1"))
index_NG<-which(NG=="1")
sampled_index<-sample(which(NG=="1"),2000) # & allgene$OG_probability< OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")

CGC_genes_index<-rep("",nrow(allgene))
CGC_genes_index[sampled_index$NG_sampled]<-"NG"
CGC_genes_index[index_TSG]<-"CGC-TSG"
CGC_genes_index[index_OG]<-"CGC-OG"
CGC_genes_index[index_TSGOG]<-"CGC-dual"
CGC_genes_index <- factor(CGC_genes_index,levels=c("CGC-dual", "CGC-OG", "CGC-TSG", "NG",""))
allgene1<-cbind(allgene,CGC_genes_index)
allgene1$CGC_genes_index <- as.factor(allgene1$CGC_genes_index)
allgene2<-allgene1[allgene1$CGC_genes_index!="",]

pdf("Raw_figures/Figure_S3E_H3K4me3_differential_methylation_scatter_for_CGC_genes.pdf", family="ArialMT", width=5, height=4.6179)#Figure 3D
p<-ggscatter(allgene2, x = "H3K4me3_peak_length", y = "Gene_body_hypermethylation_in_cancer",xlab="H3K4me3 peak length", ylab="Gene-body hypermethylation",show.legend.text=F,color = "CGC_genes_index",size = 1.6,font.label = c(9,"italic"),rug=F,label = "Core",shape=16,label.select = list(criteria = "(`x` > 3000 & `y`< 0.4) | (`y`>0.6 & `x` < 2500)"),repel=T,alpha=0.9,palette=c("darkorange2","#E41A1C","#377EB8","#CCCCCC"),alpha=0.7)+  scale_y_continuous(trans = modulus_trans(8),breaks = c(-0.4,-0.3,0,0.3,0.8),labels = c("-0.6", "-0.3", "0", "0.3", "0.8"))+ scale_x_continuous(trans = modulus_trans(1.6),breaks = c(0,2500,4000,5000,6000),labels = c("0", "2500", "4000", "5000", "8000"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
ggpar(p,legend.title = "",legend=c(0.77,0.85),font.legend = c(7, "plain", "black"))+ rremove("legend.title")
garbage <- dev.off()
