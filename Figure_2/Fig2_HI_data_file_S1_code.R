####################
##### Heatmaps #####
####################


###### DORGE predicted genes are characterized by different features.

###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/Human_TSGs_curated.txt: TSG gene list downloaded from TSGene database
###### Input: data/Human_OGs_curated.txt: OG gene list downloaded from ONGene database
###### Input: data/cancermine.txt: TSG and OG gene list downloaded from CancerMine database
###### Input: ../Figure_1/data/genebody_canyon_genes.txt: Gene-body canyon gene list
###### Output: ../Tables/Data_file_S1/Top_literature_untouched_TSG_info.txt (Data File S1): DORGE-predicted non-CGC-TSGs and related information
###### Output: ../Tables/Data_file_S1/Top_literature_untouched_OG_info.txt (Data File S1): DORGE-predicted non-CGC-OGs and related information

###### Output: Figure_2H_leftpanel.pdf,Figure_2H_rightpanel.pdf: Characterization of top 15 novel predicted genes
###### Output: Figure_2I_leftpanel.pdf,Figure_2I_rightpanel.pdf: Characterization of top 15 novel predicted genes that are undocumented in cancer driver gene databases


options(warn=-1)

### Install missing packages


installed_pkgs <- installed.packages()

pkgs <-  c("ComplexHeatmap")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
}


pkgs <-  c("circlize","dplyr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

suppressMessages(library("dplyr"))
suppressMessages(library("circlize"))
suppressMessages(library("ComplexHeatmap"))


TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005

################################## Figure 2H Top15 non-CGC TSG/OG genes ###############################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_2");

anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]
NG<-anno[,4]

prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
index_pOG<-which(prediction$OG_probability> OG_threshold & OG_CGC!="1" & TSG_CGC!="1")
index_pTSG<-which(prediction$TSG_probability>TSG_threshold & TSG_CGC!="1" & OG_CGC!="1")
all_feature <- read.table("../All_features.csv", header=T, sep=",",fill=TRUE,quote = "")
all_feature1 <- cbind(Gene=as.character(all_feature[,1]),all_feature[,-1]);
all_feature1$H3K4me3_peak_length<-as.numeric(as.character(all_feature1$H3K4me3_peak_length))
all_feature1$NonSilent_silent_ratio<-as.numeric(as.character(all_feature1$NonSilent_silent_ratio))
all_feature1$LoF_o_e_constraint<-as.numeric(as.character(all_feature1$LoF_o_e_constraint))
all_feature1$Exon_conservation_phastCons_score<-as.numeric(as.character(all_feature1$Exon_conservation_phastCons_score))
all_feature1$Frameshift_indel_fraction<-as.numeric(as.character(all_feature1$Frameshift_indel_fraction))
all_feature1$Gene_body_hypomethylation_in_cancer<-as.numeric(as.character(all_feature1$Gene_body_hypomethylation_in_cancer))
all_feature1$Missense_entropy<-as.numeric(as.character(all_feature1$Missense_entropy))
all_feature1$Missense_damaging_benign_ratio<-as.numeric(as.character(all_feature1$Missense_damaging_benign_ratio))
all_feature1$pLI_score<-as.numeric(as.character(all_feature1$pLI_score))
all_feature1$VEST_score<-as.numeric(as.character(all_feature1$VEST_score))
all_feature1$Super_enhancer_percentage<-as.numeric(as.character(all_feature1$Super_enhancer_percentage))
all_feature1$ncGERP_score<-as.numeric(as.character(all_feature1$ncGERP_score))
all_feature1$Gene_body_hypermethylation_in_cancer<-as.numeric(as.character(all_feature1$Gene_body_hypermethylation_in_cancer))
all_feature1$Gene_body_canyon_hypermethylation_in_cancer<-as.numeric(as.character(all_feature1$Gene_body_canyon_hypermethylation_in_cancer))
allgene<-all_feature1[,c("Gene","H3K4me3_peak_length","NonSilent_silent_ratio","VEST_score","Missense_damaging_benign_ratio","Exon_conservation_phastCons_score","Missense_entropy","pLI_score","Super_enhancer_percentage","ncGERP_score","Gene_body_hypermethylation_in_cancer","Gene_body_canyon_hypermethylation_in_cancer")]
allgene2<-allgene %>% 
  mutate(H3K4me3_peak_length_quantile  = ecdf(H3K4me3_peak_length)(H3K4me3_peak_length)) %>% 
  mutate(NonSilent_silent_ratio_quantile  = ecdf(NonSilent_silent_ratio)(NonSilent_silent_ratio)) %>% 
  mutate(VEST_score_quantile  = ecdf(VEST_score)(VEST_score))%>% 
  mutate(Missense_damaging_benign_ratio_quantile  = ecdf(Missense_damaging_benign_ratio)(Missense_damaging_benign_ratio))%>% 
  mutate(Exon_conservation_phastCons_score_quantile  = ecdf(Exon_conservation_phastCons_score)(Exon_conservation_phastCons_score))%>% mutate(Missense_entropy_quantile  = ecdf(Missense_entropy)(Missense_entropy))%>% 
  mutate(pLI_score_quantile  = ecdf(pLI_score)(pLI_score))%>% 
  mutate(Gene_body_canyon_hypermethylation_in_cancer_quantile  = ecdf(Gene_body_canyon_hypermethylation_in_cancer)(Gene_body_canyon_hypermethylation_in_cancer))%>% 
  mutate(Super_enhancer_percentage_quantile  = ecdf(Super_enhancer_percentage)(Super_enhancer_percentage))%>% 
  mutate(ncGERP_score_quantile  = ecdf(ncGERP_score)(ncGERP_score))%>% 
  mutate(Gene_body_hypermethylation_in_cancer_quantile  = ecdf(Gene_body_hypermethylation_in_cancer)(Gene_body_hypermethylation_in_cancer))


########## Left panel ##########
PMIDs<-c("NA","22833098","28653885","22986368","31204176","NA","NA","27121567","29342219","22538187","30655376","24570593","26993606","NA","25699089")

df_TSG<-cbind(prediction[,c("Gene","TSG_probability")],allgene2[,c("H3K4me3_peak_length_quantile","NonSilent_silent_ratio_quantile","VEST_score_quantile","Exon_conservation_phastCons_score","Super_enhancer_percentage")])
prediction_nonCGC_TSG<-df_TSG[index_pTSG,]
prediction_top_TSG<-prediction_nonCGC_TSG[order(-prediction_nonCGC_TSG[,"TSG_probability"])[c(1:15)],]#FOXP1
prediction_top_TSG<-cbind(prediction_top_TSG,PMIDs)

rownames(prediction_top_TSG)<-paste(prediction_top_TSG$Gene," PMID:",prediction_top_TSG$PMIDs,sep = "")
ha = rowAnnotation(genename = anno_text(rownames(prediction_top_TSG),just="right", location = 1.0, gp = gpar(fontsize = 7)))

pdf(file="Raw_figures/Figure_2H_leftpanel.pdf", family="ArialMT", width=3, height=5)
plot1<-ha+Heatmap(as.matrix(prediction_top_TSG[,c("H3K4me3_peak_length_quantile","NonSilent_silent_ratio_quantile","VEST_score_quantile","Exon_conservation_phastCons_score")])*100,col = colorRamp2(c(100,80,60,40,20,0),c("red", "darksalmon", "darksalmon","deepskyblue","blue","blue")),name="Quantile",show_column_dend = F,show_row_dend = F,cluster_rows=F,cluster_columns=F,show_row_names=F,show_column_names=F,heatmap_height = unit(8.5,"cm"),heatmap_width = unit(2,"cm"),heatmap_legend_param = list(title_position = "topcenter",color_bar = "continuous",legend_direction = "horizontal",legend_width = unit(1.8, "cm")),top_annotation = columnAnnotation(genename = anno_text(c("H3K4me3 peak length","Nonsilent/silent ratio","VEST score","Conservation phastCons score"),just="left", gp = gpar(fontsize = 7),location = 0,rot=45)))
draw(plot1,heatmap_legend_side = "top", padding = unit(c(10, 2, 2, 2), "mm"))
garbage <- dev.off()

########## Right panel ##########

PMIDs<-c("30367150","29629903","NA","20139090","31999936","NA","26019213","21283680","NA","27862697","28479419","31695024","29438696","28423522","22246670")

df_OG<-cbind(prediction[,c("Gene","OG_probability")],allgene2[,c("Missense_entropy_quantile","Super_enhancer_percentage_quantile","pLI_score_quantile","ncGERP_score_quantile","Gene_body_hypermethylation_in_cancer_quantile","Gene_body_canyon_hypermethylation_in_cancer_quantile")])
prediction_nonCGC_OG<-df_OG[index_pOG,]
prediction_top_OG<-prediction_nonCGC_OG[order(-prediction_nonCGC_OG[,"OG_probability"])[c(1:15)],]
prediction_top_OG<-cbind(prediction_top_OG,PMIDs)

genebody_canyon_genes <- read.table("../Figure_1/data/genebody_canyon_genes.txt", header=F)
genebody_canyon<-rep("other",nrow(prediction_top_OG))
genebody_canyon[prediction_top_OG$Gene%in%genebody_canyon_genes$V1]<-"Genebody_canyon"
prediction_top_OG$Gene_body_canyon<-genebody_canyon

rownames(prediction_top_OG)<-paste(prediction_top_OG$Gene," PMID:",prediction_top_OG$PMIDs,sep = "")

pdf(file="Raw_figures/Figure_2H_rightpanel.pdf", family="ArialMT", width=5, height=8)
ha = rowAnnotation(genename = anno_text(rownames(prediction_top_OG),just="right", location = 1, gp = gpar(fontsize = 7)))
plot1<-ha+Heatmap(as.matrix(prediction_top_OG[,c("Missense_entropy_quantile","Super_enhancer_percentage_quantile","pLI_score_quantile","ncGERP_score_quantile","Gene_body_hypermethylation_in_cancer_quantile")])*100,col = colorRamp2(c(100,80,60,40,20,0),c("red", "darksalmon", "darksalmon","deepskyblue","blue","blue")),name="DORGE-OG-score",cluster_rows=F,cluster_columns=F,show_column_dend = F,show_row_dend = F,show_row_names=F,show_column_names=F,heatmap_height = unit(8.5,"cm"),heatmap_width = unit(2,"cm"),heatmap_legend_param = list(title_position = "topcenter",color_bar = "continuous",legend_direction = "horizontal",legend_width = unit(1.8, "cm")),top_annotation = columnAnnotation(genename = anno_text(c("Missense entropy","Super enhancer percentage","pLI score","ncGERP score","Gene-body hypermethylation"),just="left", gp = gpar(fontsize = 7),location = 0,rot=45)))
draw(plot1,heatmap_legend_side = "top", padding = unit(c(10, 2, 2, 2), "mm"))

garbage <- dev.off()


################################## Figure 2I Find literature untouched TSG/OG genes ###############################

TSGene <- read.table("data/Human_TSGs_curated.txt", header=F, sep="\t")
ONGene <- read.table("data/Human_OGs_curated.txt", header=F, sep="\t")
cancermine <- read.table("data/cancermine_genelist_results.tsv", header=T, sep="\t")
cancermine_TSG<-as.character(cancermine[cancermine$tumor_suppressor_citations>0,1])
cancermine_OG<-as.character(cancermine[cancermine$oncogene_citations>0,1])


########## Left panel ##########
df_TSG<-cbind(prediction[,c("Gene","TSG_probability")],allgene2[,c("H3K4me3_peak_length_quantile","NonSilent_silent_ratio_quantile","VEST_score_quantile","Exon_conservation_phastCons_score")])
Rank_TSG<-rep(19636,19636)

TSG_qvalues<-cbind(df_TSG,anno[,4:6],Rank_TSG)
sorted = sort(TSG_qvalues[,"TSG_probability"],decreasing=T,index.return=T)
TSG_qvalues$Rank_TSG[sorted$ix]<-1:19636

TSG_sets=list(Predicted_TSG=as.character(TSG_qvalues[TSG_qvalues$TSG_probability>TSG_threshold,1]), CGC_TSG=as.character(TSG_qvalues[TSG_qvalues$TSG_all==1|TSG_qvalues$OG_all==1,1]), cancermine_TSG=cancermine_TSG,TSGene_TSG=as.character(unlist(TSGene)))

Specific_TSG <- setdiff(setdiff(setdiff(TSG_sets[[1]],TSG_sets[[2]]),TSG_sets[[3]]),TSG_sets[[4]])
TSG_selected_dat<-TSG_qvalues[which(TSG_qvalues$Gene %in% Specific_TSG),]

prediction_top_TSG<-TSG_selected_dat[order(-TSG_selected_dat[,"TSG_probability"])[c(1:15)],]
rownames(prediction_top_TSG)<-paste(prediction_top_TSG$Gene," Rank:",prediction_top_TSG$Rank_TSG,sep="")

pdf(file="Raw_figures/Figure_2I_leftpanel.pdf", family="ArialMT", width=3, height=5)
ha = rowAnnotation(genename = anno_text(rownames(prediction_top_TSG),just="right", location = 1.0, gp = gpar(fontsize = 7)))
plot1<-ha+Heatmap(as.matrix(prediction_top_TSG[,c("H3K4me3_peak_length_quantile","NonSilent_silent_ratio_quantile","VEST_score_quantile","Exon_conservation_phastCons_score")])*100,col = colorRamp2(c(100,90,80,30,0,0),c("red", "darksalmon", "darksalmon","deepskyblue","blue","blue")),name="Percentile",cluster_rows=F,cluster_columns=F,show_column_dend = F,show_row_dend = F,show_row_names=F,show_column_names=F,heatmap_height = unit(8.5,"cm"),heatmap_width = unit(2,"cm"),heatmap_legend_param = list(title_position = "topcenter",color_bar = "continuous",legend_direction = "horizontal",legend_width = unit(1.8, "cm")),top_annotation = columnAnnotation(genename = anno_text(c("H3K4me3 peak length","Nonsilent/silent ratio","VEST score","Conservation phastCons score"),just="left", gp = gpar(fontsize = 7),location = 0,rot=45)))
draw(plot1,heatmap_legend_side = "top", padding = unit(c(10, 2, 2, 2), "mm"))
garbage <- dev.off()

##########  Data file S1 Table ##########
df_TSG<-cbind(prediction[,c("Gene","TSG_probability","OG_probability")],all_feature1[,c("H3K4me3_peak_length","NonSilent_silent_ratio","VEST_score","Exon_conservation_phastCons_score")])

Rank_TSG<-rep(19636,19636)

TSG_qvalues<-cbind(df_TSG,anno[,4:6],Rank_TSG)
sorted = sort(TSG_qvalues[,2],decreasing=T,index.return=T)
TSG_qvalues$Rank_TSG[sorted$ix]<-1:19636

TSG_sets=list(Predicted_TSG=as.character(TSG_qvalues[TSG_qvalues$TSG_probability>TSG_threshold,1]), CGC_TSG=as.character(TSG_qvalues[TSG_qvalues$TSG_all==1,1]), cancermine_TSG=cancermine_TSG,TSGene_TSG=as.character(unlist(TSGene)))

Specific_TSG <- setdiff(setdiff(setdiff(TSG_sets[[1]],TSG_sets[[2]]),TSG_sets[[3]]),TSG_sets[[4]])
TSG_selected_dat<-TSG_qvalues[which(TSG_qvalues$Gene %in% Specific_TSG),]


combined<-cbind(TSG_selected_dat[,c("Gene","TSG_probability","OG_probability","H3K4me3_peak_length","NonSilent_silent_ratio","VEST_score","Exon_conservation_phastCons_score","TSG_all","OG_all","NG","Rank_TSG")])
write.table(combined,"../Tables/Data_file_S1/Top_predicted_TSG_info_without_in_databases.txt",sep="\t",quote=F,row.names=F,col.names=T)

########## Right panel ##########
df_OG<-cbind(prediction[,c("Gene","OG_probability")],allgene2[,c("Missense_entropy_quantile","Super_enhancer_percentage_quantile","pLI_score_quantile","ncGERP_score_quantile","Gene_body_hypermethylation_in_cancer_quantile","Gene_body_canyon_hypermethylation_in_cancer_quantile")])

Rank_OG<-rep(19636,19636)

OG_qvalues<-cbind(df_OG,anno[,4:6],Rank_OG)
sorted = sort(OG_qvalues[,2],decreasing=T,index.return=T)
OG_qvalues$Rank_OG[sorted$ix]<-1:19636

OG_sets=list(Predicted_OG=as.character(OG_qvalues[OG_qvalues$OG_probability> OG_threshold,1]), CGC_OG=as.character(OG_qvalues[OG_qvalues$OG_all==1|OG_qvalues$TSG_all==1,1]), cancermine_OG=cancermine_OG,ONGene_OG=as.character(unlist(ONGene)))

Specific_OG <- setdiff(setdiff(setdiff(OG_sets[[1]],OG_sets[[2]]),OG_sets[[3]]),OG_sets[[4]])
OG_selected_dat<-OG_qvalues[which(OG_qvalues$Gene %in% Specific_OG),]

prediction_top_OG<-OG_selected_dat[order(-OG_selected_dat[,"OG_probability"])[c(1:15)],]
rownames(prediction_top_OG)<-paste(prediction_top_OG$Gene," Rank:",prediction_top_OG$Rank_OG,sep="")

genebody_canyon_genes <- read.table("../Figure_1/data/genebody_canyon_genes.txt", header=F)
genebody_canyon<-rep("other",nrow(prediction_top_OG))
genebody_canyon[prediction_top_OG$Gene%in%genebody_canyon_genes$V1]<-"Genebody_canyon"
prediction_top_OG$Gene_body_canyon<-genebody_canyon
prediction_top_OG<-cbind(prediction_top_OG)


pdf(file="Raw_figures/Figure_2I_rightpanel.pdf", family="ArialMT", width=5, height=8)
ha = rowAnnotation(genename = anno_text(rownames(prediction_top_OG),just="right", location = 1, gp = gpar(fontsize = 7)))
plot1<-ha+Heatmap(as.matrix(prediction_top_OG[,4:8])*100,col = colorRamp2(c(100,90,80,30,0,0),c("red", "darksalmon", "darksalmon","deepskyblue","blue","blue")),name="Percentile",cluster_rows=F,cluster_columns=F,show_column_dend = F,show_row_dend = F,show_row_names=F,show_column_names=F,heatmap_height = unit(8.5,"cm"),heatmap_width = unit(2,"cm"),heatmap_legend_param = list(title_position = "topcenter",color_bar = "continuous",legend_direction = "horizontal",legend_width = unit(1.8, "cm")),top_annotation = columnAnnotation(genename = anno_text(c("Missense entropy","Super enhancer percentage","pLI score","ncGERP score","Gene-body hypermethylation"),just="left", gp = gpar(fontsize = 7),location = 0,rot=45)))
draw(plot1,heatmap_legend_side = "top", padding = unit(c(10, 2, 2, 2), "mm"))
garbage <- dev.off()

##########  Data file S1 Table ##########
df_OG<-cbind(prediction[,c("Gene","TSG_probability","OG_probability")],all_feature1[,c("Missense_entropy","Super_enhancer_percentage","pLI_score","ncGERP_score","Gene_body_hypermethylation_in_cancer","Gene_body_canyon_hypermethylation_in_cancer")])

Rank_OG<-rep(19636,19636)

OG_qvalues<-cbind(df_OG,anno[,4:6],Rank_OG)
sorted = sort(OG_qvalues[,3],decreasing=T,index.return=T)
OG_qvalues$Rank_OG[sorted$ix]<-1:19636

OG_sets=list(Predicted_OG=as.character(OG_qvalues[OG_qvalues$OG_probability> OG_threshold,1]), CGC_OG=as.character(OG_qvalues[OG_qvalues$OG_all==1,1]), cancermine_OG=cancermine_OG,OGene_OG=as.character(unlist(ONGene)))

Specific_OG <- setdiff(setdiff(setdiff(OG_sets[[1]],OG_sets[[2]]),OG_sets[[3]]),OG_sets[[4]])
OG_selected_dat<-OG_qvalues[which(OG_qvalues$Gene %in% Specific_OG),]

combined<-cbind(OG_selected_dat[,c("Gene","TSG_probability","OG_probability","Missense_entropy","Super_enhancer_percentage","pLI_score","ncGERP_score","Gene_body_hypermethylation_in_cancer","Gene_body_canyon_hypermethylation_in_cancer","TSG_all","OG_all","NG","Rank_OG")])
write.table(combined,"../Tables/Data_file_S1/Top_predicted_OG_info_without_in_databases.txt",sep="\t",quote=F,row.names=F,col.names=T)