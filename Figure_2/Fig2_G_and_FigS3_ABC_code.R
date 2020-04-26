###########################################
######## Heatmaps and scatterplots ########
###########################################

###### The DORGE-predicted novel TSGs and OGs are enriched in KEGG pathways. Besides, the DORGE-predicted novel TSGs and OGs using Epigenetic features only are also enriched in KEGG pathways. The difference is quantified by scatterplots.

###### Input: data/KEGG_2019_Human_table_novel_TSGs.txt: Enrichr "KEGG 2019 human" results using all DORGE-predict novel non-CGC TSGs
###### Input: data/KEGG_2019_Human_table_novel_OGs.txt: Enrichr "KEGG 2019 human" results using all DORGE-predict novel non-CGC OGs
###### Input: data/KEGG_2019_Human_table_epi_benefit_TSG.txt: Enrichr input: Gene list of TSG prediction using all features while excluding those without epigenetic features
###### Input: data/KEGG_2019_Human_table_epi_benefit_OG.txt:  Enrichr input: Gene list of OG prediction using all features while excluding those without epigenetic features
###### Input: data/KEGG_2019_Human_table_epi_benefit_TSG.txt: Enrichr "KEGG 2019" results using predicted novel TSGs based on all features while excluding those without epigenetic features
###### Input: data/KEGG_2019_Human_table_epi_benefit_OG.txt: Enrichr "KEGG 2019" results using all DORGE-predict novel OGs based on all features while excluding those without epigenetic features

###### Output: Figure_2F_KEGG_novelTSGOG.pdf: KEGG pathway enrichment using Enrichr for all DORGE-predict novel non-CGC TSG/OGs
###### Output: Figure_S3A_KEGG_benefited_by_epigenetics_novelTSGOG.pdf:  KEGG pathway enrichment using Enrichr for DORGE-predict novel TSG/OGs predicted by Epigenetics features only
###### Output: Figure_S3B_KEGG_correlation_TSG.pdf: Scatter plot comparing -log10(DORGE-TSG P-value) of KEGG pathway enrichment in Figure 2F and S4A
###### Output: Figure_S3C_KEGG_correlation_OG.pdf:  Scatter plot comparing -log10(DORGE-OG P-value) of KEGG pathway enrichment in Figure 2F and S4A

options(warn=-1)

### Install missing packages
installed_pkgs <- installed.packages()
pkgs <-  c("ComplexHeatmap")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
}

pkgs <-  c("plyr","circlize","ggpubr","ggrepel")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

suppressMessages(library(ggpubr))
suppressMessages(library(circlize))
suppressMessages(library(plyr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggrepel))

TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005

################################## Figure 2F ###############################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_2");
KEGG_TSG <- read.table("data/KEGG_2019_Human_table_novel_TSGs.txt", header=T, sep="\t")
KEGG_TSG <- KEGG_TSG[,c(1,3,4)]
names(KEGG_TSG)<-c("Term","Pvalue_T","Adjusted_T");
KEGG_OG <- read.table("data/KEGG_2019_Human_table_novel_OGs.txt", header=T, sep="\t")
KEGG_OG <- KEGG_OG[,c(1,3,4)]
names(KEGG_OG)<-c("Term","Pvalue_O","Adjusted_O");

join<-join(KEGG_TSG, KEGG_OG,type = "full",by="Term")
join$Pvalue_T[is.na(join$Pvalue_T)] <- 1
join$Adjusted_T[is.na(join$Adjusted_T)] <- 1
join$Pvalue_O[is.na(join$Pvalue_O)] <- 1
join$Adjusted_O[is.na(join$Adjusted_O)] <- 1

index_select<-which(join$Adjusted_T<0.0001| join$Adjusted_O<0.0001)
join2<-join[index_select,]

log10_join2<- cbind(join2,-log10(join2[,2:5]))
names(log10_join2)<-c("Term" ,"Pvalue_T","Adjusted_T","Pvalue_O","Adjusted_O","log10_Pvalue_T","log10_Adjusted_T","log10_Pvalue_O","log10_Adjusted_O")
index_diff<-which((log10_join2$log10_Adjusted_T-log10_join2$log10_Adjusted_O)>8| (log10_join2$log10_Adjusted_O-log10_join2$log10_Adjusted_T)>4)
index_consistent<-which(log10_join2$log10_Adjusted_T>4 & log10_join2$log10_Adjusted_O>4)
log10_join2<-log10_join2[c(which(log10_join2$Term=="Apoptosis"),union(index_diff,index_consistent)),]


dat_plot<-log10_join2[,c(1,7,9)]
dat_plot2<-dat_plot[,c(3,2)]
rownames(dat_plot2)<-dat_plot[,1]
pdf("Raw_figures/Figure_2F_KEGG_novelTSGOG.pdf", family="ArialMT", width=4, height=5.5)#Figure 2F
Heatmap(as.matrix(dat_plot2), col = colorRamp2(c(10, 5, 3, 1, 0.1), c("red","darksalmon", "deepskyblue","dodgerblue","blue")),row_names_rot = 315, row_names_side = "right",clustering_method_rows = "complete",show_row_names = T,cluster_rows=T,column_title_rot = 0,column_labels=c("OG (n=500)","TSG (n=766)"),row_dend_side = "left",row_dend_width=unit(3,"mm"),column_title_side = "top",show_column_dend = F,row_names_gp = gpar(fontsize = 7),column_names_gp = gpar(fontsize = 7),heatmap_legend_param = list(title = expression("-log"[10]~"Adjusted "~italic(P)~"-value")),heatmap_height = unit(14, "cm"),heatmap_width = unit(5.4,"cm"))
dev.off()


##### Epigenetics benefited genes ######

KEGG_TSG <- read.table("data/KEGG_2019_Human_table_epi_benefit_TSG.txt", header=T, sep="\t")
KEGG_TSG <- KEGG_TSG[,c(1,3,4)]
names(KEGG_TSG)<-c("Term","Pvalue_T","Adjusted_T");
KEGG_OG <- read.table("data/KEGG_2019_Human_table_epi_benefit_OG.txt", header=T, sep="\t")
KEGG_OG <- KEGG_OG[,c(1,3,4)]
names(KEGG_OG)<-c("Term","Pvalue_O","Adjusted_O");

join<-join(KEGG_TSG, KEGG_OG,type = "full",by="Term")
join$Pvalue_T[is.na(join$Pvalue_T)] <- 1
join$Adjusted_T[is.na(join$Adjusted_T)] <- 1
join$Pvalue_O[is.na(join$Pvalue_O)] <- 1
join$Adjusted_O[is.na(join$Adjusted_O)] <- 1

index_select<-which(join$Adjusted_T<0.01| join$Adjusted_O<0.01)
join2<-join[index_select,]

log10_join2<- cbind(join2,-log10(join2[,2:5]))
names(log10_join2)<-c("Term" ,"Pvalue_T","Adjusted_T","Pvalue_O","Adjusted_O","log10_Pvalue_T","log10_Adjusted_T","log10_Pvalue_O","log10_Adjusted_O")
#index_diff<-which((log10_join2$log10_Adjusted_T-log10_join2$log10_Adjusted_O)>2| (log10_join2$log10_Adjusted_O-log10_join2$log10_Adjusted_T)>1)
#index_consistent<-which(log10_join2$log10_Adjusted_T>2 & log10_join2$log10_Adjusted_O>2)
#log10_join2<-log10_join2[c(union(index_diff,index_consistent)),]


dat_plot_2<-log10_join2[,c(1,7,9)]
dat_plot2_2<-dat_plot_2[,c(3,2)]
rownames(dat_plot2_2)<-dat_plot_2[,1]
pdf("Raw_figures/Figure_S3A_KEGG_benefited_by_epigenetics_novelTSGOG.pdf", family="ArialMT", width=6, height=11)#Figure 2F
Heatmap(as.matrix(dat_plot2_2), col = colorRamp2(c(4, 2, 1, 0.5, 0), c("red","darksalmon", "deepskyblue","dodgerblue","blue")),row_names_rot = 315, row_names_side = "right",clustering_method_rows = "complete",show_row_names = T,cluster_rows=T,column_title_rot = 0,column_labels=c("OG (n=301)","TSG (n=564)"),row_dend_side = "left",row_dend_width=unit(3,"mm"),column_title_side = "top",show_column_dend = F,row_names_gp = gpar(fontsize = 7),column_names_gp = gpar(fontsize = 7),heatmap_legend_param = list(title = expression("-log"[10]~"Adjusted "~italic(P)~"-value")),heatmap_height = unit(10, "cm"),heatmap_width = unit(5.055,"cm"))
dev.off()


##### Figure S3B: Correlation of KEGG enrichment P-values between novel TSGs predicted on all features and without epigenetic features ######

KEGG_TSG <- read.table("data/KEGG_2019_Human_table_novel_TSGs.txt", header=T, sep="\t")
KEGG_TSG <- KEGG_TSG[,c(1,3,4)]
names(KEGG_TSG)<-c("Term","Pvalue_T","Adjusted_T");
KEGG_OG <- read.table("data/KEGG_2019_Human_table_novel_OGs.txt", header=T, sep="\t")
KEGG_OG <- KEGG_OG[,c(1,3,4)]
names(KEGG_OG)<-c("Term","Pvalue_O","Adjusted_O");

join<-join(KEGG_TSG, KEGG_OG,type = "full",by="Term")
join$Pvalue_T[is.na(join$Pvalue_T)] <- 1
join$Adjusted_T[is.na(join$Adjusted_T)] <- 1
join$Pvalue_O[is.na(join$Pvalue_O)] <- 1
join$Adjusted_O[is.na(join$Adjusted_O)] <- 1

index_select<-which(join$Adjusted_T<0.01| join$Adjusted_O<0.01)
join2<-join[index_select,]

log10_join2<- cbind(join2,-log10(join2[,2:5]))

dat_plot<-log10_join2[,c(1,7,9)]
dat_plot2<-dat_plot[,c(3,2)]
rownames(dat_plot2)<-dat_plot[,1]


KEGG_TSG <- read.table("data/KEGG_2019_Human_table_epi_benefit_TSG.txt", header=T, sep="\t")
KEGG_TSG <- KEGG_TSG[,c(1,3,4)]
names(KEGG_TSG)<-c("Term","Pvalue_T","Adjusted_T");
KEGG_OG <- read.table("data/KEGG_2019_Human_table_epi_benefit_OG.txt", header=T, sep="\t")
KEGG_OG <- KEGG_OG[,c(1,3,4)]
names(KEGG_OG)<-c("Term","Pvalue_O","Adjusted_O");

join<-join(KEGG_TSG, KEGG_OG,type = "full",by="Term")
join$Pvalue_T[is.na(join$Pvalue_T)] <- 1
join$Adjusted_T[is.na(join$Adjusted_T)] <- 1
join$Pvalue_O[is.na(join$Pvalue_O)] <- 1
join$Adjusted_O[is.na(join$Adjusted_O)] <- 1

index_select<-which(join$Adjusted_T<0.05| join$Adjusted_O<0.05)
join2<-join[index_select,]

log10_join2<- cbind(join2,-log10(join2[,2:5]))
names(log10_join2)<-c("Term" ,"Pvalue_T","Adjusted_T","Pvalue_O","Adjusted_O","log10_Pvalue_T","log10_Adjusted_T","log10_Pvalue_O","log10_Adjusted_O")



dat_plot_2<-log10_join2[,c(1,7,9)]
dat_plot2_2<-dat_plot_2[,c(3,2)]
rownames(dat_plot2_2)<-dat_plot_2[,1]

names(dat_plot)<-c("Term","log10_Pvalue_T","log10_Pvalue_O");
names(dat_plot_2)<-c("Term","log10_Pvalue_T2","log10_Pvalue_O2");
combined<-join(dat_plot,dat_plot_2,type = "left",by="Term")

pdf("Raw_figures/Figure_S3B_KEGG_correlation_TSG.pdf", family="ArialMT", width=4, height=4)

suppressMessages(print(ggscatter(combined, x = "log10_Pvalue_T", y = "log10_Pvalue_T2",add = "reg.line",add.params = list(color = "purple", fill = "purple"),xlab="-log10 P-value of KEGG enrichment\nfor DORGE-predicted novel TSGs",ylab="-log P-value of KEGG\nenrichment for Epigenetic feature\nbenefited DORGE-predicted novel TSGs") +geom_text_repel(data = combined[combined$log10_Pvalue_T2 >2.5,],aes(label=Term,x=log10_Pvalue_T, y=log10_Pvalue_T2), family="ArialMT",size= 3,color="purple",segment.size= 0.3,force = 15,segment.color= "grey10") + stat_cor(method ="spearman",cor.coef.name="rho",digits=2,label.sep="\n", label.x = 3, label.y = 3.5))) # Add correlation coefficient
dev.off()

##### Figure S3C: Correlation of KEGG enrichment P-values between novel OGs predicted on all features and without epigenetic features ######
pdf("Raw_figures/Figure_S3C_KEGG_correlation_OG.pdf", family="ArialMT", width=4, height=4)
suppressMessages(print(ggscatter(combined, x = "log10_Pvalue_O", y = "log10_Pvalue_O2",add = "reg.line",xlab="-log P-value of KEGG enrichment\nfor DORGE-predicted novel OGs",ylab="-log10 P-value of KEGG\nenrichment for Epigenetic feature\nbenefited DORGE-predicted novel OGs",add.params = list(color = "red", fill = "red"))+geom_text_repel(data = combined[combined$log10_Pvalue_O2 >3.5 & combined$log10_Pvalue_O2 <6,],aes(label=Term,x=log10_Pvalue_O, y=log10_Pvalue_O2), family="ArialMT",size= 3,color="red",segment.size= 0.3,force = 15,segment.color= "grey10")+ylim(0,8)+stat_cor(method ="spearman",cor.coef.name="rho",digits=2,label.sep="\n", label.x = 10, label.y = 5))) # Add correlation coefficient
dev.off()
