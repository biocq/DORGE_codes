#####################################
##### Scatterplots and boxplots #####
#####################################

###### The scatterplots of Missense entropy and gene-body differential methylation (and also gene-body differential methylation at gene-body canyon genes) are generated, highlighted by DORGE-OG scores.

###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: sampled_NGs.txt: Sampled NGs for reproducing purpose
###### Input: data/genebody_canyon_genes.txt: Gene-body canyon genes (Compiled from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1492-3/figures/5)


###### Output: Figure_S3F_Missense_entropy_genebody_hyper_methyl_scatter.pdf: Scatterplot of Missense entropy and gene-body differential methylation
###### Output: Figure_S3G_Missense_entropy_and_diff_methyl_scatter.pdf: Scatterplot of Missense entropy and gene-body differential methylation at canyon genes
###### Output: Figure_S3H_Missense_entropy_boxplot.pdf: Missense entropy comparison for gene-body canyon genes, NGs and CGC-OGs

options(warn=-1)

### Install missing packages


installed_pkgs <- installed.packages()

pkgs <-  c("ggpubr","scales","dplyr","cowplot")


if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}


suppressMessages(library("dplyr"))
suppressMessages(library("ggpubr"))
suppressMessages(library("scales"))
suppressMessages(library("cowplot"))

TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005

############### Figure S3F: Scatter plot of gene-body differential methylation and Missense entropy #############

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
NG<-anno[,"NG"]
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","PTEN","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","FOXP1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1"
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
mis_entropy <- read.table("data/Missense_entropy_accurate_version.txt", header=T, sep="\t",fill=TRUE,quote = "")
all_feature<-all_feature[,c("Gene","Gene_body_hypermethylation_in_cancer")]
all_feature<-cbind(all_feature,mis_entropy[,c("Missense_entropy")])
colnames(all_feature)<-c("Gene","Gene_body_hypermethylation_in_cancer","Missense_entropy")
lowdiff<-which(all_feature$Gene_body_hypermethylation_in_cancer< -0.3)
all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]<-(all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]+0.3)/3-0.3
index_highentropy<-which(all_feature$Missense_entropy>2)
all_feature$Missense_entropy <- as.numeric(all_feature$Missense_entropy)
all_feature$Missense_entropy[index_highentropy]<-(all_feature$Missense_entropy[index_highentropy]-2)/2+2


allgene<-cbind(allgene,all_feature[,c("Gene_body_hypermethylation_in_cancer","Missense_entropy")])
allgene$Missense_entropy[is.na(allgene$Missense_entropy)] <- 0
allgene$Gene_body_hypermethylation_in_cancer[is.na(allgene$Gene_body_hypermethylation_in_cancer)] <- 0


index<-rep("",nrow(allgene))
index_E<-which(allgene$OG_probability> OG_threshold)
#allgene$Missense_entropy>0.15 & allgene$Gene_body_hypermethylation_in_cancer<0.6 & 
index_H<-which(allgene$OG_probability> OG_threshold)
#allgene$Gene_body_hypermethylation_in_cancer>0 & allgene$Missense_entropy<1 & 

index[index_E]<-"Missense_entropy"
index[index_H]<-"GB_hyper"
index_eh<-c(index_E,index_H)
allgene<-cbind(allgene,index)
allgene$index <- as.character(allgene$index)
sampled_index<-sample(which(NG=="1" & allgene[,"index"]==""),2000)# & allgene$OG_probability< OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")
allgene[sampled_index$NG_sampled,"index"] <-"Other"

logP_OG<- as.numeric(allgene$OG_probability)
combined_logP<-rep(0,nrow(allgene))
combined_logP[index_H]<-logP_OG[index_H]
combined_logP[index_E]<-logP_OG[index_E]
allgene1<-cbind(allgene,combined_logP)
allgene2<-allgene1[allgene1$index!="",]


pdf("Raw_figures/Figure_S3F_Missense_entropy_genebody_hyper_methyl_scatter.pdf", family="ArialMT", width=5, height=5)#Figure 3C
p<-ggscatter(allgene2, x = "Gene_body_hypermethylation_in_cancer", y = "Missense_entropy",xlab="Gene-body differential methylation", ylab="Missense entropy",show.legend.text=F,color = "combined_logP", size = 1.2, font.label = c(9, "italic"),label = "Core",shape=16, repel = T,label.select = list(criteria = "`combined_logP` > 0.9 & (`x` > 0.75 | `y` > 1.5)"),alpha=0.9)+ gradient_color(rev(c("#FF0000","#DA6969","#DA6969","#DA6969","#D37E7E","#D37E7E","#D37E7E","#D37E7E","#CC9393","#CC9393","#CC9393","#CC9393","#CC9393","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#BEBEBE","#BEBEBE","#BEBEBE")))+  scale_y_continuous(trans = modulus_trans(-1),breaks = c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,2,4),labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.6", "0.8", "1", "2", "6"))+ scale_x_continuous(trans = modulus_trans(6),breaks = c(-0.4,-0.3,0,0.3,0.8),labels = c("-0.6", "-0.3", "0", "0.3", "0.8"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "transparent",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))

ggpar(p,legend.title = "",legend=c(0.89,0.8),legend.background = element_rect(color = "white",fill = "", size = 1, linetype = "solid"),ylim=c(0,6),font.legend = c(7, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=0.78, y=2.3, size=2, label="DORGE-OG score",color="red",angle = 90)

garbage <- dev.off()

####### Figure S3G: Scatter plot of gene-body differential methylation vs Missense entropy at gene-body canyon genes #######




#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_3");

anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
NG<-anno[,"NG"]
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core")]
allgene$Core <-rep("",19636);
allgene$Gene <- as.character(allgene$Gene)
allgene[which(allgene$Gene %in% c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","PTEN","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","FOXP1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1")),"Core"]<-"1"
allgene[which(allgene$Core=="1"),"Core"]<-as.character(allgene$Gene[which(allgene$Core=="1")])
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")

mis_entropy <- read.table("data/Missense_entropy_accurate_version.txt", header=T, sep="\t",fill=TRUE,quote = "")

all_feature<-all_feature[,c("Gene","Gene_body_hypermethylation_in_cancer")]
all_feature<-cbind(all_feature,mis_entropy[,c("Missense_entropy")])
lowdiff<-which(all_feature$Gene_body_hypermethylation_in_cancer< -0.3)
all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]<-(all_feature$Gene_body_hypermethylation_in_cancer[lowdiff]+0.3)/3-0.3

colnames(all_feature)<-c("Gene","Gene_body_hypermethylation_in_cancer","Missense_entropy")
index_highentropy<-which(all_feature$Missense_entropy>2)
all_feature$Missense_entropy <- as.numeric(all_feature$Missense_entropy)
all_feature$Missense_entropy[index_highentropy]<-(all_feature$Missense_entropy[index_highentropy]-2)/2+2


allgene<-cbind(allgene,all_feature[,c("Gene_body_hypermethylation_in_cancer","Missense_entropy")])
allgene$Missense_entropy[is.na(allgene$Missense_entropy)] <- 0
allgene$Gene_body_hypermethylation_in_cancer[is.na(allgene$Gene_body_hypermethylation_in_cancer)] <- 0

genebody_canyon_genes<-read.table("../Figure_1/data/genebody_canyon_genes.txt", header=F)
genebody_canyon_gene_list<-as.character(genebody_canyon_genes$V1)
index_canyon<-which(allgene$Gene_body_hypermethylation_in_cancer>-10 & all_feature$Gene %in% genebody_canyon_gene_list)
index<-rep("",nrow(allgene))
index[index_canyon]<-"GB_hyper_canyon"
allgene<-cbind(allgene,index)
allgene$index <- as.character(allgene$index)

sampled_index<-sample(which(NG=="1"),2000) # & allgene[,"index"]=="" & allgene$OG_probability< OG_threshold & allgene$TSG_probability<TSG_threshold
sampled_index <- read.table("sampled_NGs.txt", header=T, sep="\t",fill=TRUE,quote = "")
#allgene[sampled_index$NG_sampled,"index"] <-"Other"

logP_OG<- as.numeric(allgene$OG_probability)
combined_logP<-rep(0,nrow(allgene))
combined_logP[index_canyon]<-logP_OG[index_canyon]
allgene1<-cbind(allgene,combined_logP)
allgene2<-allgene1[allgene1$index!="" & allgene1$Gene_body_hypermethylation_in_cancer>-10.3,]
allgene3<-allgene2;
add<-allgene2[1,]
add["Missense_entropy"]<-4.5
add["Gene"]<-"Not_used"
add["Gene_body_hypermethylation_in_cancer"]<- -100
allgene3<-rbind(allgene3,add)

pdf("Raw_figures/Figure_S3G_Missense_entropy_and_diff_methyl_scatter.pdf", family="ArialMT", width=5, height=5)#Figure 3E
p<-ggscatter(allgene3, x = "Gene_body_hypermethylation_in_cancer", y = "Missense_entropy",xlab="Differential methylation at gene-body canyon", ylab="Missense entropy",show.legend.text=F,color = "combined_logP", size = 1.2, font.label = c(9, "italic"),label = "Core",shape=16, repel = T,label.select = list(criteria = "`combined_logP` > 0"),alpha=0.9)+ gradient_color(rev(c("#FF0000","#DA6969","#DA6969","#DA6969","#D37E7E","#D37E7E","#D37E7E","#D37E7E","#CC9393","#CC9393","#CC9393","#CC9393","#CC9393","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#C5A8A8","#BEBEBE","#BEBEBE","#BEBEBE")))+  scale_y_continuous(trans = modulus_trans(-1),breaks = c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,2,4),labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.6", "0.8", "1", "2", "6"))+ scale_x_continuous(trans = modulus_trans(10),breaks = c(-0.3,0,0.3,0.8),labels = c("-0.3","0", "0.3", "0.8"))+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "transparent",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))

ggpar(p,legend.title = "",legend=c(0.3,0.8),legend.background = element_rect(color = "white",fill = "", size = 1, linetype = "solid"),xlim=c(-0.4,0.8),font.legend = c(7, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=0.42, y=2, size=2, label="DORGE-OG score",color="red",angle = 90)+ annotate(geom="text", x=0.572, y=1.47, size=2.5, label="0.00",color="black")

garbage <- dev.off()

############ Figure S3H: Boxplot of missense entropy of gene-body canyon genes, CGC-OGs and NGs #############

TSG_CGC<-anno[,"TSG_all"]
OG_CGC<-anno[,"OG_all"]
NG<-anno[,"NG"]

index_CGC_OG<-which(OG_CGC=="1")
index_NG<-which(NG=="1")
canyon_genes<-rep("",nrow(allgene))
canyon_genes[index_canyon]<-"Canyon"
allgene<-cbind(allgene,canyon_genes)
index_NG<-which(NG=="1")
gene_cat<-rep("",nrow(allgene))
gene_cat[index_NG]<-"NG"
gene_cat[index_CGC_OG]<-"OG"
gene_cat[index_canyon]<-"Canyon"
allgene<-cbind(allgene,gene_cat)     
allgene2<-allgene[allgene$gene_cat!="",]

my_comparisons <- list( c("NG", "Canyon"),c("OG", "Canyon"))

pdf("Raw_figures/Figure_S3H_Missense_entropy_boxplot.pdf", family="ArialMT", width=2, height=5.6927)#Figure 3FG
ggboxplot(allgene2, x = "gene_cat", y = "Missense_entropy", xlab="", ylab="Missense entropy", color = "gene_cat", fill = "gene_cat",outlier.size = 0.2,width=0.7, palette =c("orange","grey","red"),alpha = 0.7)+ theme_bw()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(hjust = 0.5, colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+  scale_y_continuous(trans = modulus_trans(-1),breaks = c(0,0.1,0.2,0.3,0.4,0.6,0.8,1,2,4),labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.6", "0.8", "1", "2", "6"))+ guides(fill=FALSE)+ border()+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "g"),comparisons = my_comparisons,size=2.5)+ rremove("legend")
garbage <- dev.off()
