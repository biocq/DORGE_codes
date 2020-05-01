########################
##### Scatterplots #####
########################

###### The scatterplots of H3K4me3 peak length and Missense entropy (and also differential gene-body methylation) are generated, highlighted by DORGE-predicted genes with high H3K4me3 peak length or high Missense entropy.

###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: sampled_NGs.txt: Sampled NGs for reproducing purpose
###### Input: Missense_entropy_accurate_version.txt: Missense entropy feature with more accurate values

###### Output: Figure_3C_H3K4me3_Missense_Entropy_scatter.pdf: Scatterplot of H3K4me3 peak length and Missense entropy
###### Output: Figure_3E_H3K4me3_diff_methyl_scatter.pdf Scatterplot of H3K4me3 peak length and differential gene-body methylation

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
##################### Figure 3C: Scatter plot of Broad H3K4me3 peak length and Missense entropy #####################

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

mis_entropy <- read.table("data/Missense_entropy_accurate_version.txt", header=T, sep="\t",fill=TRUE,quote = "")
anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
NG<-anno[,"NG"]
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","PTEN","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","FOXP1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1"
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
allgene<-cbind(allgene,all_feature[,c("H3K4me3_peak_length")],mis_entropy[,c("Missense_entropy")])
colnames(allgene)<-c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","Core","H3K4me3_peak_length","Missense_entropy")
allgene$Missense_entropy[is.na(allgene$Missense_entropy)] <- 0
allgene$H3K4me3_peak_length[is.na(allgene$H3K4me3_peak_length)] <- 0
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>5000)
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-5000)/3+5000

index_highentropy<-which(allgene$Missense_entropy>2)
allgene$Missense_entropy <- as.numeric(allgene$Missense_entropy)
allgene$Missense_entropy[index_highentropy]<-(allgene$Missense_entropy[index_highentropy]-2)/2+2

index<-rep("",nrow(allgene))
index_E<-which(allgene$Missense_entropy>0 & allgene$H3K4me3_peak_length<4000 & allgene$OG_probability>OG_threshold)
index_H<-which(allgene$H3K4me3_peak_length>0 & allgene$Missense_entropy<0.6 & allgene$TSG_probability>TSG_threshold)
index[index_E]<-"Missense_entropy"
index[index_H]<-"Broad_H3K4me3"
index_eh<-c(index_E,index_H)
allgene<-cbind(allgene,index)
allgene$index <- as.character(allgene$index)
sampled_index<-sample(which(NG=="1"),2000) # & allgene[,"index"]=="" & allgene$OG_probability<OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")
allgene[sampled_index$NG_sampled,"index"] <-"Other"
minuslogP_TSG<- as.numeric(allgene$TSG_probability)
logP_OG<- -as.numeric(allgene$OG_probability)
combined_logP<-rep(0,nrow(allgene))
combined_logP[index_H]<-minuslogP_TSG[index_H]
combined_logP[index_E]<-logP_OG[index_E]
allgene1<-cbind(allgene,combined_logP)
allgene1$combined_logP <- as.numeric(allgene1$combined_logP)
allgene1$index <- as.factor(allgene1$index)
allgene2<-allgene1[allgene1$index!="",]


pdf("Raw_figures/Figure_3C_H3K4me3_Missense_Entropy_scatter.pdf", family="ArialMT", width=5, height=5)#Figure 3C
p<-ggscatter(allgene2, x = "H3K4me3_peak_length", y = "Missense_entropy",xlab="H3K4me3 peak length", ylab="Missense entropy",show.legend.text=F,color = "combined_logP", size = 1.6, font.label = c(9, "italic"),label = "Core",shape=16, repel = T,label.select = list(criteria = "`combined_logP` > 0.99| `combined_logP` < -0.999999"),alpha=0.9)+ gradient_color(c("#FF0000","#DA6969","#DA6969","#D37E7E","#D37E7E","#D37E7E","#CC9393","#CC9393","#CC9393","#CC9393","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#BEBEBE","#BEBEBE","#BEBEBE","#A8A8C5","#A8A8C5","#9393CC","#9393CC","#9393CC","#9393CC","#9393CC","#7E7ED3","#7E7ED3","#7E7ED3","#7E7ED3","#6969DA","#6969DA","#6969DA","#5454E2","#5454E2","#0000FF"))+  scale_y_continuous(trans = modulus_trans(-1.5),breaks = c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,2,4),labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.6", "0.8", "1", "2", "6"))+ scale_x_continuous(trans = modulus_trans(1.6),breaks = c(0,2500,4000,5000,6000),labels = c("0", "2500", "4000", "5000", "8000"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "transparent",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
ggpar(p,legend.title = "",legend=c(0.78,0.8),legend.background = element_rect(color = "transparent",fill = "transparent", size = 1, linetype = "solid"),ylim=c(0,6),font.legend = c(7, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=5000, y=3.5, size=2, label="DORGE-TSG score",color="blue",angle = 90)+ annotate(geom="text", x=5000, y=1, size=2, label="DORGE-OG score",color="red",angle = 90)+ annotate(geom="text", x=5540, y=2.76, size=2.5, label="1.0",color="black")+geom_rect(aes(xmin=5500, xmax=5560, ymin=1.1, ymax=1.5), fill='grey4',alpha=1)
garbage <- dev.off()

########### Figure 3E: Scatter plot of Broad H3K4me3 peak length and gene-body differential methylation  ############

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_3");

all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","H3K4me3_peak_length","Gene_body_hypermethylation_in_cancer")]

lowdiff<-which(all_feature$Gene_body_hypermethylation_in_cancer< -0.3)
all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]<-(all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]+0.3)/3-0.3
anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
NG<-anno[,"NG"]
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","PTEN","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1"#"FOXP1",
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])
allgene<-cbind(allgene,all_feature[,c("H3K4me3_peak_length","Gene_body_hypermethylation_in_cancer")])

allgene$Gene_body_hypermethylation_in_cancer[is.na(allgene$Gene_body_hypermethylation_in_cancer)] <- 0
allgene$H3K4me3_peak_length[is.na(allgene$H3K4me3_peak_length)] <- 0
index_BroadK4me3<-which(allgene$H3K4me3_peak_length>5000)
allgene$H3K4me3_peak_length <- as.numeric(allgene$H3K4me3_peak_length)
allgene$H3K4me3_peak_length[index_BroadK4me3]<-(allgene$H3K4me3_peak_length[index_BroadK4me3]-5000)/3+5000

index<-rep("",nrow(allgene))
index_H<-which(allgene$Gene_body_hypermethylation_in_cancer>0 & allgene$H3K4me3_peak_length<4000 & allgene$OG_probability> OG_threshold)
index_B<-which(allgene$H3K4me3_peak_length>0 & allgene$Gene_body_hypermethylation_in_cancer<0.6 & allgene$TSG_probability>TSG_threshold)
index[index_B]<-"Broad_H3K4me3"
index[index_H]<-"GB_hyper"
index_eh<-c(index_H,index_B)
allgene<-cbind(allgene,index)
allgene$index <- as.character(allgene$index)
sampled_index<-sample(which(NG=="1"),2000) # & allgene[,"index"]=="" & allgene$OG_probability <OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")
allgene[sampled_index$NG_sampled,"index"] <-"Other"

minuslogP_TSG<- as.numeric(allgene$TSG_probability)
logP_OG<- -as.numeric(allgene$OG_probability)
combined_logP<-rep(0,nrow(allgene))
combined_logP[index_B]<-minuslogP_TSG[index_B]
combined_logP[index_H]<-logP_OG[index_H]
allgene1<-cbind(allgene,combined_logP)
allgene1$combined_logP <- as.numeric(allgene1$combined_logP)
allgene1$index <- as.factor(allgene1$index)
allgene2<-allgene1[allgene1$index!="",]

pdf("Raw_figures/Figure_3E_H3K4me3_diff_methyl_scatter.pdf", family="ArialMT", width=5, height=5)#Figure S4D
p<-ggscatter(allgene2, x = "H3K4me3_peak_length", y = "Gene_body_hypermethylation_in_cancer",xlab="H3K4me3 peak length", ylab="Gene-body hypermethylation",show.legend.text=F,color = "combined_logP", size = 1.6, font.label = c(9, "italic"),label = "Core",shape=16, repel = T,label.select = list(criteria = "(`combined_logP` > 0.99 & `x` > 4000)| (`combined_logP` < -0.9999999 & `y` > 0.6)"),alpha=0.9)+ gradient_color(c("#FF0000","#DA6969","#DA6969","#D37E7E","#D37E7E","#D37E7E","#CC9393","#CC9393","#CC9393","#CC9393","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#BEBEBE","#BEBEBE","#BEBEBE","#A8A8C5","#A8A8C5","#9393CC","#9393CC","#9393CC","#9393CC","#9393CC","#7E7ED3","#7E7ED3","#7E7ED3","#7E7ED3","#6969DA","#6969DA","#6969DA","#5454E2","#5454E2","#0000FF"))+  scale_y_continuous(trans = modulus_trans(5),breaks = c(-0.4,-0.3,0,0.3,0.8),labels = c("-0.6", "-0.3", "0", "0.3", "0.8"))+ scale_x_continuous(trans = modulus_trans(1.6),breaks = c(0,2500,4000,5000,6000),labels = c("0", "2500", "4000", "5000", "8000"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "transparent",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))

ggpar(p,legend.title = "",legend=c(0.9,0.8),legend.background = element_rect(color = "transparent",fill = "transparent", size = 1, linetype = "solid"),ylim=c(-0.41,0.8),font.legend = c(7, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=5500, y=0.76, size=2, label="DORGE-TSG score",color="blue",angle = 90)+ annotate(geom="text", x=5500, y=0.65, size=2, label="DORGE-OG score",color="red",angle = 90)
garbage <- dev.off()
