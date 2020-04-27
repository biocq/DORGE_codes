####################
##### Heatmaps #####
####################

###### Heatmap of fold changes in expression for DORGE-predicted TSG/OGs that are upregulated or downregulated by drug perturbations.

###### Input data are available at Evaluation_processing folder DORGE_pipeline project https://github.com/biocq/DORGE_pipeline
###### Input: data/DGB/DGB_creeds_downregulated_simplified.txt: Processed data including fold changes in expression for DORGE-predicted OGs (excluding CGC dual genes and DORGE-predicted dual genes) that are downregulated by drug perturbations from the CRowd Extracted Expression of Differential Signatures (CREEDS) collection in the Drug Gene Budger (DGB) web portal. Only most signficant drugs and genes are shown due to the limitation of the space.
###### Input: data/DGB/DGB_creeds_upregulated_simplified.txt: Processed data including fold changes in expression for DORGE-predicted TSG (excluding CGC dual genes and DORGE-predicted dual genes) that are upregulated by drug perturbations from the CRowd Extracted Expression of Differential Signatures (CREEDS) collection in the Drug Gene Budger (DGB) web portal. Only most signficant drugs and genes are shown due to the limitation of the space.
###### Input: data/DGB/DGB_CMap_downregulated_simplified.txt: Processed data including fold changes in expression for DORGE-predicted OGs (excluding CGC dual genes and DORGE-predicted dual genes) that are downregulated by drug perturbations from the Connectivity Map (CMap) collection in the Drug Gene Budger (DGB) web portal. Only most signficant drugs and genes are shown due to the limitation of the space.
###### Input: data/DGB/DGB_CMap_upregulated_simplified.txt: Processed data including fold changes in expression for DORGE-predicted TSGs (excluding CGC dual genes and DORGE-predicted dual genes) that are upregulated by drug perturbations from the Connectivity Map (CMap) collection in the Drug Gene Budger (DGB) web portal. Only most signficant drugs and genes are shown due to the limitation of the space.
###### Output: Figure_S6_creeds_heatmap.pdf: Heatmap of fold changes in expression for DORGE-predicted TSG/OGs that are up-regulated or down-regulated by drug perturbations.
###### Output: Figure_S7_CMap_heatmap.pdf: Heatmap of fold changes in expression for DORGE-predicted TSG/OGs that are up- or down-regulated by drug perturbations.

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("data.table","circlize","plyr","Hmisc","reshape2")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

pkgs <-  c("ComplexHeatmap")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
}

suppressMessages(library("plyr"))
suppressMessages(library("data.table"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("Hmisc"))
suppressMessages(library("reshape2"))

############################# Figure S6 Heatmap of compound-gene heatmap based on CREEDS data #########################

nonzero<- function(x) {
   length(x[x!=0])
}
plus<- function(x) {
   length(x[x>0])
}
minus<- function(x) {
   length(x[x<0])
}

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_5");
creeds_dataset <- read.table("data/DGB/DGB_creeds_downregulated_simplified.txt", header=T, sep="\t",quote = "")
creeds_dataset2 <- read.table("data/DGB/DGB_creeds_upregulated_simplified.txt", header=T, sep="\t",quote = "")
creeds_dat<-rbind(creeds_dataset,creeds_dataset2)
creeds_dat_2fc<-creeds_dat[(creeds_dat$Fold_change>1) | (creeds_dat$Fold_change< -1),]#log2 fold change indeed
gene_count<-count(creeds_dat_2fc, "Gene")
recurrent_genes<-as.character(gene_count[gene_count$freq>0,1]) # Only genes with compounds of frequency >0 are shown
drug_count<-count(creeds_dat_2fc[which(creeds_dat_2fc$Gene %in% recurrent_genes),], "Compound")
recurrent_drugs<-as.character(drug_count[drug_count$freq>0,1]) # Only compounds with genes of frequency >0 are shown

creeds_dat_2fc_bak<-creeds_dat_2fc[which(creeds_dat_2fc$Compound %in% recurrent_drugs),]
creeds_dat_2fc_bak<-creeds_dat_2fc_bak[which(creeds_dat_2fc_bak$Gene %in% recurrent_genes),]
dat_wide<-creeds_dat_2fc_bak[,c(1,2,5,6)]
cast.data<-reshape2::dcast(dat_wide, Gene+Gene_type~Compound, fun.aggregate = mean,value.var="Fold_change")

Gene_type<-as.character(cast.data$Gene_type)
dat_plot<-cast.data[,-2]
rownames(dat_plot)<-cast.data[,1]
dat_plot<-as.matrix(dat_plot[,-1])
dat_plot[is.nan(dat_plot)] <- 0

dat_plot<-dat_plot[order(Gene_type),]
count<-as.numeric(apply(dat_plot,1,FUN=nonzero))#Gene
Gene_type<-Gene_type[order(Gene_type)]
rowselect<-which(count>6)
Gene_type<-Gene_type[rowselect]
dat_plot<-dat_plot[rowselect,]

count<-as.numeric(apply(dat_plot,2,FUN=nonzero))#Drug
colselect<-which(count>2)
colnames(dat_plot)<-capitalize(colnames(dat_plot))

compound_names <-c("3,3',4,4'-tetrachlorobiphenyl"="Compounds with cell-line evidence","4-hydroxynonenal"="Compounds with cell-line evidence","Actinomycind"="Anti-cancer and chemotherapy drugs","Adenosinetriphosphate"="Compounds without anti-cancer evidence","Androstanolone"="Compounds without anti-cancer evidence","Antimony"="Compounds without anti-cancer evidence","Apratoxina"="Compounds with cell-line evidence","Arachidonicacid"="Compounds with cell-line evidence","Azacitidine"="Anti-cancer and chemotherapy drugs","Bexarotene"="Anti-cancer and chemotherapy drugs","Bisphenola"="Compounds without anti-cancer evidence","Bpde"="Compounds without anti-cancer evidence","Calcitriol"="Anti-cancer and chemotherapy drugs","Chromium"="Compounds without anti-cancer evidence","Cisplatin"="Anti-cancer and chemotherapy drugs","Clinafloxacin"="Anti-cancer and chemotherapy drugs","Curcumin"="Compounds with cell-line evidence","Cycloheximide"="Compounds with cell-line evidence","Cytarabine"="Anti-cancer and chemotherapy drugs","Dasatinib"="Anti-cancer and chemotherapy drugs","Dexamethasone"="Anti-cancer and chemotherapy drugs","Diclofenac"="Anti-cancer and chemotherapy drugs","Dmnq"="Compounds with cell-line evidence","Doxorubicin"="Anti-cancer and chemotherapy drugs","Doxycycline"="Anti-cancer and chemotherapy drugs","Estradiol"="Compounds without anti-cancer evidence","Etanercept"="Compounds without anti-cancer evidence","Gefitinib"="Anti-cancer and chemotherapy drugs","Hydrocortisone"="Anti-cancer and chemotherapy drugs","Imatinib"="Anti-cancer and chemotherapy drugs","Interferongamma-1b"="Anti-cancer and chemotherapy drugs","Lapatinib"="Anti-cancer and chemotherapy drugs","Lipopolysaccharide"="Compounds with cell-line evidence","Luteolin"="Anti-cancer and chemotherapy drugs","Lysophosphatidicacid"="Compounds without anti-cancer evidence","Mercury"="Compounds without anti-cancer evidence","Metformin"="Compounds with cell-line evidence","Methotrexate"="Anti-cancer and chemotherapy drugs","Mln4924"="Compounds with cell-line evidence","Naturalalphainterferon"="Anti-cancer and chemotherapy drugs","Neocarzinostatin"="Anti-cancer and chemotherapy drugs","Nickel"="Compounds without anti-cancer evidence","Nilotinib"="Anti-cancer and chemotherapy drugs","Nitricoxide"="Compounds with cell-line evidence","Pd173074"="Compounds with cell-line evidence","Phorbol12-myristate13-acetate"="Compounds with cell-line evidence","PMA"="Compounds with cell-line evidence","Plx4032"="Anti-cancer and chemotherapy drugs","Plx4720"="Compounds with cell-line evidence","Propofol"="Compounds with cell-line evidence","Puromycin,2xec50,3d"="Compounds with cell-line evidence","Puromycin,ec50,1d"="Compounds with cell-line evidence","Puromycin,ec50,3d"="Compounds with cell-line evidence","Puromycin,ec50,5d"="Compounds with cell-line evidence","R1881"="Compounds with cell-line evidence","Resveratrol"="Compounds with cell-line evidence","Ribavirin"="Anti-cancer and chemotherapy drugs","Sapphyrinpci-2050"="Anti-cancer and chemotherapy drugs","Sevoflurane"="Compounds with cell-line evidence","Sulforaphane"="Anti-cancer and chemotherapy drugs","Thapsigargin"="Compounds with cell-line evidence","Torcetrapib"="Compounds with cell-line evidence","Tretinoin"="Anti-cancer and chemotherapy drugs","Trovafloxacin"="Compounds with cell-line evidence","Ubiquinol"="Compounds without anti-cancer evidence","Vanadiumpentoxide"="Compounds without anti-cancer evidence","Vemurafenib"="Anti-cancer and chemotherapy drugs","Vx"="Compounds without anti-cancer evidence","Y15"="Compounds with cell-line evidence")
compound_names <- c(compound_names,"15-deltaprostaglandinJ2"="Compounds with cell-line evidence","alvespimycin"="Anti-cancer and chemotherapy drugs","cephaeline"="Compounds with cell-line evidence","cloperastine"="Compounds without anti-cancer evidence","diethylstilbestrol"="Compounds without anti-cancer evidence","digoxigenin"="Compounds without anti-cancer evidence","geldanamycin"="Anti-cancer and chemotherapy drugs","helveticoside"="Compounds with cell-line evidence","ionomycin"="Compounds with cell-line evidence","lanatosideC"="Compounds with cell-line evidence","lycorine"="Compounds with cell-line evidence","phenoxybenzamine"="Compounds with cell-line evidence","rosiglitazone"="Compounds with cell-line evidence","tanespimycin"="Anti-cancer and chemotherapy drugs","thioridazine"="Compounds with cell-line evidence","tretinoin"="Anti-cancer and chemotherapy drugs","trichostatinA"="Compounds with cell-line evidence","trifluoperazine"="Compounds with cell-line evidence","troglitazone"="Compounds with cell-line evidence")
compound_names <- c(compound_names,"Phorbol12-myristate13-acetate(pma)"="Compounds with cell-line evidence","Promegestone"="Anti-cancer and chemotherapy drugs","Aplidin"="Anti-cancer and chemotherapy drugs","Bicalutamide"="Anti-cancer and chemotherapy drugs","Quercetin"="Compounds with cell-line evidence","Carboplatin(36h)"="Anti-cancer and chemotherapy drugs","Tibolone"="Compounds without anti-cancer evidence","Chlorpyrifos"="Compounds without anti-cancer evidence","Tamoxifen"="Anti-cancer and chemotherapy drugs","Dmnq(2,3-dimethoxy-1,4-naphthoquinone)"="Compounds without anti-cancer evidence","Sapphyrinpci-2050(1.25&icirc;&frac14;m)"="Compounds with cell-line evidence","Sapphyrinpci-2050(2.5&icirc;&frac14;m)"="Compounds with cell-line evidence","Pioglitazone"="Compounds with cell-line evidence","Fulvestrant"="Anti-cancer and chemotherapy drugs","Genistein"="Anti-cancer and chemotherapy drugs","Valproicacid"="Anti-cancer and chemotherapy drugs","beta-escin"="Compounds with cell-line evidence","Gossypol"="Compounds with cell-line evidence","Monorden"="Compounds with cell-line evidence","Prochlorperazine"="Compounds with cell-line evidence")
names(compound_names)<-capitalize(names(compound_names))



pdf("Raw_figures/Figure_S6_creeds_heatmap.pdf", family="ArialMT", width=11.32, height=8.1)
ha = rowAnnotation(Gene_type = Gene_type, col = list(Gene_type = c("CGC-OG" = "#E41A1C","Novel DORGE-OG" = "#4DAF4A", "CGC-TSG" = "#377EB8", "Novel DORGE-TSG" = "#984EA3")), show_annotation_name = F)
ha2 = HeatmapAnnotation("Gene number" = anno_barplot(cbind(as.numeric(apply(dat_plot[,colselect],2,FUN=plus)),as.numeric(apply(dat_plot[,colselect],2,FUN=minus))),gp = gpar(fill = c("blue","red"), col = c("blue","red"),border="white")))
ha3 = columnAnnotation(Compound_annotation = compound_names[colnames(dat_plot[,colselect])], col = list(Compound_annotation = c("Anti-cancer and chemotherapy drugs" = "#CCCCFF","Compounds with cell-line evidence" = "#FFCC99", "Compounds without anti-cancer evidence" = "#00CCCC")), show_annotation_name = F,simple_anno_size = unit(0.15, "cm"))

Heatmap(dat_plot[,colselect], col = colorRamp2(c(-3,-1,0,1,3), c("red","pink","white","lightblue","blue")),show_row_names = T,name="Log2 fold change", row_names_gp = gpar(fontsize = 7),column_names_gp = gpar(fontsize = 7),cluster_rows=T,cluster_columns=T,show_column_dend=F, show_row_dend=F,heatmap_legend_param = list(title_position = "leftcenter-rot",color_bar = "continuous"), right_annotation =ha, top_annotation =ha2, bottom_annotation =ha3)
garbage <- dev.off()

############################# Figure S7 Heatmap of compound-gene heatmap based on CMap data #########################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_5");
CMap_dataset <- read.table("data/DGB/DGB_CMap_downregulated_simplified.txt", header=T, sep="\t",quote = "")
CMap_dataset2 <- read.table("data/DGB/DGB_CMap_upregulated_simplified.txt", header=T, sep="\t",quote = "")
CMap_dat<-rbind(CMap_dataset,CMap_dataset2)
CMap_dat_2fc<-CMap_dat[(CMap_dat$Fold_change>1) | (CMap_dat$Fold_change< -1),]
gene_count<-count(CMap_dat_2fc, "Gene")
recurrent_genes<-as.character(gene_count[gene_count$freq>0,1])
drug_count<-count(CMap_dat_2fc[which(CMap_dat_2fc$Gene %in% recurrent_genes),], "Compound")
recurrent_drugs<-as.character(drug_count[drug_count$freq>0,1])

CMap_dat_2fc_bak<-CMap_dat_2fc[which(CMap_dat_2fc$Compound %in% recurrent_drugs),]
CMap_dat_2fc_bak<-CMap_dat_2fc_bak[which(CMap_dat_2fc$Gene %in% recurrent_genes),]
dat_wide<-CMap_dat_2fc_bak[,c(1,2,5,6)]
cast.data<-reshape2::dcast(dat_wide, Gene+Gene_type~Compound, fun.aggregate = mean,value.var="Fold_change")

Gene_type<-as.character(cast.data$Gene_type)
dat_plot<-cast.data[,-2]
rownames(dat_plot)<-cast.data[,1]
dat_plot<-as.matrix(dat_plot[,-1])
dat_plot[is.nan(dat_plot)] <- 0

dat_plot<-dat_plot[order(Gene_type),]
count<-as.numeric(apply(dat_plot,1,FUN=nonzero))#Gene
Gene_type<-Gene_type[order(Gene_type)]
rowselect<-which(count>2)
Gene_type<-Gene_type[rowselect]
dat_plot<-dat_plot[rowselect,]

count<-as.numeric(apply(dat_plot,2,FUN=nonzero))#Drug
#dat_plot<-dat_plot[,rev(order(as.numeric(count)))]
colselect<-which(count>1)
colnames(dat_plot)<-capitalize(colnames(dat_plot))


pdf("Raw_figures/Figure_S7_CMap_heatmap.pdf", family="ArialMT", width=7.935, height=7.26)

ha = rowAnnotation(Gene_type = Gene_type, col = list(Gene_type = c("CGC-OG" = "#E41A1C","Novel DORGE-OG" = "#4DAF4A", "CGC-TSG" = "#377EB8", "Novel DORGE-TSG" = "#984EA3")), show_annotation_name = F)
ha2 = HeatmapAnnotation("Gene number" = anno_barplot(cbind(as.numeric(apply(dat_plot[,colselect],2,FUN=plus)),as.numeric(apply(dat_plot[,colselect],2,FUN=minus))),gp = gpar(fill = c("blue","red"), col = c("blue","red"),border="white")))
ha3 = columnAnnotation(Compound_annotation = compound_names[colnames(dat_plot[,colselect])], col = list(Compound_annotation = c("Anti-cancer and chemotherapy drugs" = "#CCCCFF","Compounds with cell-line evidence" = "#FFCC99", "Compounds without anti-cancer evidence" = "#00CCCC")), show_annotation_name = F,simple_anno_size = unit(0.35, "cm"))
#fill according to colnames(dat_plot[,colselect])
colnames(dat_plot)<-capitalize(colnames(dat_plot))
Heatmap(dat_plot[,colselect], col = colorRamp2(c(-3,-1,0,1,3), c("red","pink","white","lightblue","blue")),show_row_names = T,name="Log2 fold change", row_names_gp = gpar(fontsize = 7),column_names_gp = gpar(fontsize = 7),cluster_rows=T,cluster_columns=T,show_column_dend=F, show_row_dend=F,heatmap_legend_param = list(title_position = "leftcenter-rot",color_bar = "continuous"), right_annotation =ha, top_annotation =ha2, bottom_annotation =ha3)
garbage <- dev.off()