####################
##### QQ-plots #####
####################

###### DORGE predictions are shown by QQ-plots with specific gene categoried highlighted.

###### Input: data/TSG_prediction-lr_enet_wo_Epigenetics_FPR_0.01.csv: DORGE-TSG predicted genes using features without epigenetic features
###### Input: data/OG_prediction-lr_enet_wo_Epigenetics_FPR_0.01.csv: DORGE-OG prediction genes using features without epigenetic features
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Output: Figure_3A_qqplot_TSG_with_allfeature_specific_highlighted.pdf: QQplot for DORGE-TSG prediction
###### Output: Figure_3B_qqplot_OG_with_allfeature_specific_highlighted.pdf: QQplot for DORGE-OG prediction

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggpubr","scales","ggrepel")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("ggpubr"))
suppressMessages(library("scales"))
suppressMessages(library("ggrepel"))

TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

##################### Figure 3A #####################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_3");
score_tsg <- read.table("data/TSG_prediction-lr_enet_wo_Epigenetics_FPR_0.01.csv", header=T, sep=",")
Predicted_genes_woepi_features<-as.character(score_tsg[,1])
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")

index_predicted_novel_TSGs_allfeatures<-which(prediction$TSG_probability>TSG_threshold)
Novel_genes_allfeatures<-as.character(prediction[index_predicted_novel_TSGs_allfeatures,1])
Specific_novel_genes_allfeatures<-setdiff(Novel_genes_allfeatures,Predicted_genes_woepi_features)
CGC_TSGs<-as.character(prediction[prediction$TSG_all=="1",1])
Specific_novel_nonCGC_genes_allfeatures<-setdiff(Specific_novel_genes_allfeatures,CGC_TSGs)
write.table(Specific_novel_nonCGC_genes_allfeatures,"data/Novel_DORGE_predicted_TSGs_benefitted_by_Epi_features.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Specific_novel_genes_allfeatures,"data/DORGE_predicted_TSGs_benefitted_by_Epi_features.txt",sep="\t",quote=F,row.names=F,col.names=F)

index_TSG_predicted_CGC_TSGs_all_feature<-which(prediction$TSG_probability>TSG_threshold & prediction$TSG_all=="1")
all_genes<-as.character(prediction[,1])
Predicted_CGC_TSGs_all_features<-all_genes[index_TSG_predicted_CGC_TSGs_all_feature]
index_predicted_woepi_genes_in_all_genes<-which(all_genes %in% Predicted_genes_woepi_features)
Predicted_CGC_TSGs_nonepi_features_only<-all_genes[intersect(index_predicted_woepi_genes_in_all_genes,index_TSG_predicted_CGC_TSGs_all_feature)]
all_feature_specific_genes<-setdiff(Predicted_CGC_TSGs_all_features,Predicted_CGC_TSGs_nonepi_features_only)

index_unpredicted_CGC_TSGs<-which(prediction$TSG_probability<TSG_threshold & prediction$TSG_all=="1")
unpredicted_CGC_TSGs<-all_genes[index_unpredicted_CGC_TSGs]

sorted = sort(prediction$TSG_probability,decreasing=F,index.return=T)
sorted_genes<-as.character(all_genes[sorted$ix])
index_highlighted_epi_specific_genes<-which(sorted_genes %in% all_feature_specific_genes)
index_highlighted_nonepi_specific_genes<-which(sorted_genes %in% Predicted_CGC_TSGs_nonepi_features_only)
index_highlighted_unpredicted_genes<-which(sorted_genes %in% unpredicted_CGC_TSGs)
highlighted_genes<-sorted_genes[c(which(sorted_genes %in% c(all_feature_specific_genes,Predicted_CGC_TSGs_nonepi_features_only,unpredicted_CGC_TSGs)))]

nonhighlighted<-as.character(sorted_genes[which(! sorted_genes %in% highlighted_genes)])

other_selected<- which(sorted_genes %in% nonhighlighted)

o <- sorted$x
e <-  1:length(prediction$TSG_probability)/19636
e[1]<-0
epi_specific<-rep("",19636)
epi_specific[index_highlighted_epi_specific_genes]<-"Epigenetic_specific"
nonepi_specific<-rep("",19636)
nonepi_specific[index_highlighted_nonepi_specific_genes]<-"woEpigenetic_specific"
unpredicted_CGC_TSGs<-rep("",19636)
unpredicted_CGC_TSGs[index_highlighted_unpredicted_genes]<-"Unpredicted CGC-TSGs"
colors<-rep("",19636)
colors[index_highlighted_epi_specific_genes]<-"Epigenetic_specific"
colors[index_highlighted_nonepi_specific_genes]<-"woEpigenetic_specific"
colors[index_highlighted_unpredicted_genes]<-"Unpredicted CGC-TSGs"

core_genes<-c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","PTEN","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","FOXP1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1");

core_gene_vec<-rep("",19636)
core_gene_vec[sorted_genes %in% core_genes]<-"Core_CGC"

dat <- data.frame("Observed" = o, "Expected" = e,"Epi_specific_genes" = epi_specific,"nonEpi_specific_genes" = nonepi_specific,"Unpredicted_CGC" = unpredicted_CGC_TSGs,"Core" = core_gene_vec,"Gene" = sorted_genes,"Color" = colors)



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

dat1<-dat
dat<-dat[dat$Color!="",]
label<-which(dat[,3]!="" & dat[,1]>TSG_threshold & dat[,1]<0.9999)
label_chosen<-label[seq(1, length(label), 5)]

pdf("Raw_figures/Figure_3A_qqplot_TSG_with_allfeature_specific_highlighted.pdf", family="ArialMT", width=3.9, height=2.7)
dat$Color<-factor(dat$Color,levels=c("Epigenetic_specific","woEpigenetic_specific","Unpredicted CGC-TSGs"))
palette_color <- c("Epigenetic_specific"='red', "woEpigenetic_specific"='blue', "Unpredicted CGC-TSGs"='darkgreen')

ggplot() + geom_path(data=dat1, aes(x=Expected, y=Observed))+ geom_point(data=dat[dat[,"Color"]!="",], aes(Expected, Observed, colour=Color),position=position_jitter(w = 0.1, h = 0.1, seed = 2),size=0.3)+ scale_colour_manual(values = palette_color,breaks = c("Epigenetic_specific","woEpigenetic_specific","Unpredicted CGC-TSGs"),labels = c(paste("CGC-TSGs (n=",length(index_highlighted_epi_specific_genes),") benefit from\nepigenetic features only",sep=""), paste("CGC-TSGs (n=",length(index_highlighted_nonepi_specific_genes),") can be predicted\n by non-epigenetic features",sep=""), paste("CGC-TSGs (n=",length(index_highlighted_unpredicted_genes),") not predicted\nby DORGE-TSG",sep="")))+geom_text_repel(data = dat[label_chosen,],aes(label=Gene,x=Expected, y=Observed), family="ArialMT",size= 3,color="red",segment.size= 0.3,force = 1,segment.color= "grey10")+ geom_hline(yintercept=TSG_threshold, linetype="dashed", color = "grey10")+annotate("text", label = paste("FPR = 0.01, TSG-score = ",round(TSG_threshold,4),sep=""), x = 0.45, y =  TSG_threshold+0.02, size = 2, colour = "black")+scale_x_continuous(trans = modulus_trans(5),breaks = c(0,TSG_threshold,0.8,0.9,0.95,1),labels = c("0", "Threshold","0.8","0.9", "0.95", "1"))+scale_y_continuous(trans = modulus_trans(5),breaks = c(0,TSG_threshold,0.8,0.9,0.95,1),labels = c("0", "Threshold","0.8","0.9", "0.95", "1"))+ xlab("Expected TSG-score") + ylab("Observed TSG-score")+ guides(color = guide_legend(override.aes=list(fill=NA))) + theme_bw() + theme(legend.title = element_blank(),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_line(size = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.4,0.8),legend.background = element_rect(fill=alpha('transparent', 0.0)))

garbage <- dev.off()


##################### Figure 3B #####################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_3");
score_OG <- read.table("data/OG_prediction-lr_enet_wo_Epigenetics_FPR_0.01.csv", header=T, sep=",")
Predicted_genes_woepi_features<-as.character(score_OG[,1])
prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")

index_predicted_novel_OGs_allfeatures<-which(prediction$OG_probability>OG_threshold)
Novel_genes_allfeatures<-as.character(prediction[index_predicted_novel_OGs_allfeatures,1])
Specific_novel_genes_allfeatures<-setdiff(Novel_genes_allfeatures,Predicted_genes_woepi_features)
CGC_OGs<-as.character(prediction[prediction$OG_all=="1",1])
Specific_novel_nonCGC_genes_allfeatures<-setdiff(Specific_novel_genes_allfeatures,CGC_OGs)
write.table(Specific_novel_nonCGC_genes_allfeatures,"data/Novel_DORGE_predicted_OGs_benefitted_by_Epi_features.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Specific_novel_genes_allfeatures,"data/DORGE_predicted_OGs_benefitted_by_Epi_features.txt",sep="\t",quote=F,row.names=F,col.names=F)

index_OG_predicted_CGC_OGs_all_feature<-which(prediction$OG_probability>OG_threshold & prediction$OG_all=="1")
all_genes<-as.character(prediction[,1])
Predicted_CGC_OGs_all_features<-all_genes[index_OG_predicted_CGC_OGs_all_feature]
index_predicted_woepi_genes_in_all_genes<-which(all_genes %in% Predicted_genes_woepi_features)
Predicted_CGC_OGs_nonepi_features_only<-all_genes[intersect(index_predicted_woepi_genes_in_all_genes,index_OG_predicted_CGC_OGs_all_feature)]
all_feature_specific_genes<-setdiff(Predicted_CGC_OGs_all_features,Predicted_CGC_OGs_nonepi_features_only)

index_unpredicted_CGC_OGs<-which(prediction$OG_probability<OG_threshold & prediction$OG_all=="1")
unpredicted_CGC_OGs<-all_genes[index_unpredicted_CGC_OGs]

sorted = sort(prediction$OG_probability,decreasing=F,index.return=T)
sorted_genes<-as.character(all_genes[sorted$ix])
index_highlighted_epi_specific_genes<-which(sorted_genes %in% all_feature_specific_genes)
index_highlighted_nonepi_specific_genes<-which(sorted_genes %in% Predicted_CGC_OGs_nonepi_features_only)
index_highlighted_unpredicted_genes<-which(sorted_genes %in% unpredicted_CGC_OGs)
highlighted_genes<-sorted_genes[c(which(sorted_genes %in% c(all_feature_specific_genes,Predicted_CGC_OGs_nonepi_features_only,unpredicted_CGC_OGs)))]

nonhighlighted<-as.character(sorted_genes[which(! sorted_genes %in% highlighted_genes)])

other_selected<- which(sorted_genes %in% nonhighlighted)

o <- sorted$x
e <-  1:length(prediction$OG_probability)/19636
e[1]<-0
epi_specific<-rep("",19636)
epi_specific[index_highlighted_epi_specific_genes]<-"Epigenetic_specific"
nonepi_specific<-rep("",19636)
nonepi_specific[index_highlighted_nonepi_specific_genes]<-"woEpigenetic_specific"
unpredicted_CGC_OGs<-rep("",19636)
unpredicted_CGC_OGs[index_highlighted_unpredicted_genes]<-"Unpredicted CGC-OGs"
colors<-rep("",19636)
colors[index_highlighted_epi_specific_genes]<-"Epigenetic_specific"
colors[index_highlighted_nonepi_specific_genes]<-"woEpigenetic_specific"
colors[index_highlighted_unpredicted_genes]<-"Unpredicted CGC-OGs"

core_genes<-c("MYC","RET","TP53","BRAF","ERBB2","KRAS","ABL1","KIT","PTEN","MYCN","WT1","RB1","NRAS","ALK","CDKN2A","BCR","MET","APC","FOS","BRCA1","BCL2","CCND1","VHL","YAP1","AKT1","MYB","NOTCH1","RUNX3","KLF4","EWSR1","PIK3CA","KMT2A","CDH1","STAT3","HRAS","EZH2","JUN","ROS1","SRC","SHH","GLI1","BCL6","CTNNB1","PML","STK11","RASSF1","HGF","SMAD4","NF1","SOX2","MEN1","MARK2","LMO2","CDKN1A","PDGFRA","RUNX1","MITF","MTDH","FOXM1","MLLT3","MTOR","ERG","NPM1","MDM2","TAZ","EGF","CADM1","BRCA2","PDGFB","CD274","VEGFA","TERT","AR","CSF1R","MYCL","FHIT","KLF6","SYK","PTPN11","PTTG1","NF2","NTRK1","WWOX","FGFR1","TCL1A","ARID1A","DLC1","FOXP1","CDK16","PTGS2","RARA","SLC3A2","BAP1","FBXW7","TAL1","CDK8","ETS1","PTGDR","TNFAIP3","TWSG1");

core_gene_vec<-rep("",19636)
core_gene_vec[sorted_genes %in% core_genes]<-"Core_CGC"

dat <- data.frame("Observed" = o, "Expected" = e,"Epi_specific_genes" = epi_specific,"nonEpi_specific_genes" = nonepi_specific,"Unpredicted_CGC" = unpredicted_CGC_OGs,"Core" = core_gene_vec,"Gene" = sorted_genes,"Color" = colors)

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

dat1<-dat
dat<-dat[dat$Color!="",]
label<-which(dat[,3]!="" & dat[,1]>OG_threshold & dat[,1]<0.9999)
label_chosen<-label[seq(1, length(label), 5)]

pdf("Raw_figures/Figure_3B_qqplot_OG_with_allfeature_specific_highlighted.pdf", family="ArialMT", width=3.9, height=2.7)
dat$Color<-factor(dat$Color,levels=c("Epigenetic_specific","woEpigenetic_specific","Unpredicted CGC-OGs"))
palette_color <- c("Epigenetic_specific"='red', "woEpigenetic_specific"='blue', "Unpredicted CGC-OGs"='darkgreen')

ggplot() + geom_path(data=dat1, aes(x=Expected, y=Observed))+ geom_point(data=dat[dat[,"Color"]!="",], aes(Expected, Observed, colour=Color),position=position_jitter(w = 0.1, h = 0.1, seed = 2),size=0.3)+ scale_colour_manual(values = palette_color,breaks = c("Epigenetic_specific","woEpigenetic_specific","Unpredicted CGC-OGs"),labels = c(paste("CGC-OGs (n=",length(index_highlighted_epi_specific_genes),") benefit from\nepigenetic features only",sep=""), paste("CGC-OGs (n=",length(index_highlighted_nonepi_specific_genes),") can be predicted\n by non-epigenetic features",sep=""), paste("CGC-OGs (n=",length(index_highlighted_unpredicted_genes),") not predicted\nby DORGE-OG",sep="")))+geom_text_repel(data = dat[label_chosen,],aes(label=Gene,x=Expected, y=Observed), family="ArialMT",size= 3,color="red",segment.size= 0.3,force = 1,segment.color= "grey10")+ geom_hline(yintercept=OG_threshold, linetype="dashed", color = "grey10")+annotate("text", label = paste("FPR = 0.01, OG-score = ",round(OG_threshold,4),sep=""), x = 0.42, y =  OG_threshold+0.02, size = 2, colour = "black")+scale_x_continuous(trans = modulus_trans(5),breaks = c(0,OG_threshold,0.8,0.9,0.95,1),labels = c("0", "Threshold","0.8","0.9", "0.95", "1"))+scale_y_continuous(trans = modulus_trans(5),breaks = c(0,OG_threshold,0.8,0.9,0.95,1),labels = c("0", "Threshold","0.8","0.9", "0.95", "1"))+ xlab("Expected OG-score") + ylab("Observed OG-score")+ guides(color = guide_legend(override.aes=list(fill=NA))) + theme_bw() + theme(legend.title = element_blank(),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_line(size = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=c(0.4,0.8),legend.background = element_rect(fill=alpha('transparent', 0.0)))

garbage <- dev.off()