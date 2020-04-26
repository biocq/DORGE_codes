########################
##### Scatterplots #####
########################

###### The comparison of H3K4me3 peak length and Missense entropy for DORGE-predicted non-CGC genes and also all CGC genes are compared.

options(warn=-1)

###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: sampled_NGs.txt: Sampled NGs for reproducing purpose

###### Output: Figure_3D_H3K4me3_Missense_Entropy_scatter_for_DORGE_novel_genes.pdf
###### Output: Figure_S3D_H3K4me3_Missense_Entropy_scatter_for_CGC_genes.pdf


### Install missing packages


installed_pkgs <- installed.packages()

pkgs <-  c("ggpubr","scales","dplyr","circlize")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

suppressMessages(library("circlize"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggpubr"))
suppressMessages(library("scales"))

TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005

######### Figure 3D: Scatter plot of Broad H3K4me3 peak length and Missense entropy for DORGE-predicted novel genes #########

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

anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
TSG_CGC<-anno[,"TSG_all"]
OG_CGC<-anno[,"OG_all"]
NG<-anno[,"NG"]

allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","TSG_all","OG_all","NG")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1" #,"PTEN","FOXP1"
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])

mis_entropy <- read.table("data/Missense_entropy_accurate_version.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-cbind(allgene,all_feature[,c("H3K4me3_peak_length")],mis_entropy[,c("Missense_entropy")])
colnames(allgene)<-c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","TSG_all","OG_all","NG","Core","H3K4me3_peak_length","Missense_entropy")
allgene$Missense_entropy[is.na(allgene$Missense_entropy)] <- 0
allgene$H3K4me3_peak_length[is.na(allgene$H3K4me3_peak_length)] <- 0
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>5000)
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-5000)/3+5000

index_highentropy<-which(allgene$Missense_entropy>2)
allgene$Missense_entropy <- as.numeric(allgene$Missense_entropy)
allgene$Missense_entropy[index_highentropy]<-(allgene$Missense_entropy[index_highentropy]-2)/2+2

index<-rep("",nrow(allgene))
index_TSG<-which((TSG_CGC=="1" & OG_CGC!="1"))
index_OG<-which((OG_CGC=="1" & TSG_CGC!="1"))
index_TSGOG<-which((OG_CGC=="1" & TSG_CGC=="1"))
index_NG<-which(NG=="1")
index_pTSGOG<-which((allgene$TSG_probability>TSG_threshold & OG_CGC=="1" & TSG_CGC!="1")|(allgene$OG_probability> OG_threshold & TSG_CGC=="1" & OG_CGC!="1")| (allgene$TSG_probability>TSG_threshold & allgene$OG_probability> OG_threshold & OG_CGC!="1" & TSG_CGC!="1"))
index_pOG<-c(which(allgene$OG_probability> OG_threshold & allgene$TSG_probability< TSG_threshold & OG_CGC!="1"))
index_pTSG<-which(allgene$TSG_probability> TSG_threshold & allgene$OG_probability< OG_threshold & TSG_CGC!="1")

sampled_index<-sample(which(NG=="1"),2000)# & allgene$OG_probability< OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")

DORGE_novel_index<-rep("",nrow(allgene))
DORGE_novel_index[sampled_index$NG_sampled]<-"NG"
DORGE_novel_index[index_pTSG]<-"Novel DORGE-TSG"
DORGE_novel_index[index_pOG]<-"Novel DORGE-OG"
DORGE_novel_index[index_pTSGOG]<-"Novel DORGE-dual"
DORGE_novel_index <- factor(DORGE_novel_index,levels=c("Novel DORGE-dual", "Novel DORGE-OG", "Novel DORGE-TSG", "NG",""))
allgene1<-cbind(allgene,DORGE_novel_index)

allgene2<-allgene1[allgene1$DORGE_novel_index!="",]


pdf("Raw_figures/Figure_3D_H3K4me3_Missense_Entropy_scatter_for_DORGE_novel_genes.pdf", family="ArialMT", width=5, height=4.6179)#Figure 3D
p<-ggscatter(allgene2, x = "H3K4me3_peak_length", y = "Missense_entropy",xlab="H3K4me3 peak length", ylab="Missense entropy",show.legend.text=F,color = "DORGE_novel_index",size = 1.6,font.label = c(9,"italic"),rug=F,label = "Core",shape=16,label.select = list(criteria = "`x` > 4000 | `y`>2"),repel=T,alpha=0.9,palette=c("darkorange2","#4DAF4A","#984EA3","#CCCCCC"),alpha=0.7)+  scale_y_continuous(trans = modulus_trans(-1.5),breaks = c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,2,4),labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.6", "0.8", "1", "2", "6"))+ scale_x_continuous(trans = modulus_trans(1.6),breaks = c(0,2500,4000,5000,6000),labels = c("0", "2500", "4000", "5000", "8000"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
ggpar(p,legend.title = "",legend=c(0.85,0.87),font.legend = c(7, "plain", "black"))+ rremove("legend.title")
garbage <- dev.off()

######### Figure S3D: Scatter plot of Broad H3K4me3 peak length and Missense entropy for CGC genes #########
anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
TSG_CGC<-anno[,"TSG_all"]
OG_CGC<-anno[,"OG_all"]
NG<-anno[,"NG"]

allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","TSG_all","OG_all","NG")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1" #,"PTEN","FOXP1"
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])

mis_entropy <- read.table("data/Missense_entropy_accurate_version.txt", header=T, sep="\t",fill=TRUE,quote = "")

all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-cbind(allgene,all_feature[,c("H3K4me3_peak_length")],mis_entropy[,c("Missense_entropy")])

colnames(allgene)<-c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","TSG_all","OG_all","NG","Core","H3K4me3_peak_length","Missense_entropy")

allgene$Missense_entropy[is.na(allgene$Missense_entropy)] <- 0
allgene$H3K4me3_peak_length[is.na(allgene$H3K4me3_peak_length)] <- 0
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>5000)
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-5000)/3+5000

index_highentropy<-which(allgene$Missense_entropy>2)
allgene$Missense_entropy <- as.numeric(allgene$Missense_entropy)
allgene$Missense_entropy[index_highentropy]<-(allgene$Missense_entropy[index_highentropy]-2)/2+2

index<-rep("",nrow(allgene))
index_TSG<-which((TSG_CGC=="1" & OG_CGC!="1"))
index_OG<-which((OG_CGC=="1" & TSG_CGC!="1"))
index_TSGOG<-which((OG_CGC=="1" & TSG_CGC=="1"))
index_NG<-which(NG=="1")
sampled_index<-sample(which(NG=="1"),2000)# & allgene$OG_probability<OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")

CGC_genes_index<-rep("",nrow(allgene))
CGC_genes_index[sampled_index$NG_sampled]<-"NG"
CGC_genes_index[index_TSG]<-"CGC-TSG"
CGC_genes_index[index_OG]<-"CGC-OG"
CGC_genes_index[index_TSGOG]<-"CGC-dual"
CGC_genes_index <- factor(CGC_genes_index,levels=c("CGC-dual", "CGC-OG", "CGC-TSG", "NG",""))
allgene1<-cbind(allgene,CGC_genes_index)



allgene2<-allgene1[allgene1$CGC_genes_index!="",]
pdf("Raw_figures/Figure_S3D_H3K4me3_Missense_Entropy_scatter_for_CGC_genes.pdf", family="ArialMT", width=5, height=4.6179)#Figure S4E
p<-ggscatter(allgene2, x = "H3K4me3_peak_length", y = "Missense_entropy",xlab="H3K4me3 peak length", ylab="Missense entropy",show.legend.text=F,color = "CGC_genes_index",size = 1.6,font.label = c(9,"italic"),rug=F,label = "Core",shape=16,label.select = list(criteria = "(`x` > 3000 & `y`< 0.6) | (`y`>0.8 & `y`< 1.8 & `x` < 2500)"),repel=T,alpha=0.9,palette=c("darkorange2","#E41A1C","#377EB8","#CCCCCC"),alpha=0.7)+  scale_y_continuous(trans = modulus_trans(-1.5),breaks = c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,2,4),labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.6", "0.8", "1", "2", "6"))+ scale_x_continuous(trans = modulus_trans(1.6),breaks = c(0,2500,4000,5000,6000),labels = c("0", "2500", "4000", "5000", "8000"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))

ggpar(p,legend.title = "",legend=c(0.85,0.87),font.legend = c(7, "plain", "black"))+ rremove("legend.title")

garbage <- dev.off()
