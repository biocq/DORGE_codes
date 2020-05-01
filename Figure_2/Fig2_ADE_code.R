#########################
##### Venn diagrams #####
#########################

###### The DORGE-predicted TSGs and OGs are compared with the different databases (CGC, TSGene, ONGene and CancerMine).

###### Input: data/Human_TSGs_curated.txt: TSG gene list downloaded from TSGene database
###### Input: data/Human_OGs_curated.txt: OG gene list downloaded from ONGene database
###### Input: data/cancermine.txt: TSG and OG gene list downloaded from CancerMine database
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: ../Figure_1/data/genebody_canyon_genes.txt: Gene-body canyon gene list
###### Input: ../Gene_set_new.txt: Gene annotation file
###### Output: Figure_2A_DORGE_compared_with_CGC_venn.pdf: Venn diagram comparing DORGE prediction with CGC gene annotation
###### Output: Figure_2D_DORGE_TSG_cancermine_TSGene_venn.pdf: Venn diagram comparing DORGE TSG prediction with CancerMine and TSGene gene lists
###### Output: Figure_2E_DORGE_OG_cancermine_ONGene_venn.pdf: Venn diagram comparing DORGE OG prediction with ONGene gene list


options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("VennDiagram")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library(VennDiagram))


TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005
######  Data preparation  ######

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_2");
TSGene <- read.table("data/Human_TSGs_curated.txt", header=F, sep="\t")
ONGene <- read.table("data/Human_OGs_curated.txt", header=F, sep="\t")

anno <- read.table("../Gene_set_new.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(anno)<-c("Gene","TSG_core","OG_core","NG","TSG_all","OG_all");
TSG_CGC<-anno[,5]
OG_CGC<-anno[,6]

prediction <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")

index_pOG<-which(prediction$OG_probability> OG_threshold)
index_pTSG<-which(prediction$TSG_probability>TSG_threshold)
index_TSG<-which(TSG_CGC=="1")
index_OG<-which(OG_CGC=="1")


new_predicted_ogs<-data.frame("Gene"=anno[index_pOG,1])
known_og_genes<-data.frame("Gene"=anno[index_OG,1])
new_predicted_tsgs<-data.frame("Gene"=anno[index_pTSG,1])
known_tsg_genes<-data.frame("Gene"=anno[index_TSG,1])


###### Figure 2A: DORGE-TSGs, DORGE-OGs, CGC-TSGs and CGC-OGs ######
TSGOG_sets=list("DORGE-TSGs"=as.character(new_predicted_tsgs[,1]), "CGC-OGs"=as.character(known_og_genes[,1]), "CGC-TSGs"=as.character(known_tsg_genes[,1]), "DORGE-OGs"=as.character(new_predicted_ogs[,1]))
temp <- venn.diagram(TSGOG_sets, fill=c("#984EA3","#E41A1C","#377EB8","#4DAF4A"), alpha=c(0.7,0.7,0.7,0.7), cex=1, cat.fontface=4,filename=NULL)

pdf(file="Raw_figures/Figure_2A_DORGE_compared_with_CGC_venn.pdf", family="ArialMT", width=4, height=4)
    grid.draw(temp)
dev.off()

###### Figure 2D: DORGE-TSGs, Cancermine and TSGene ######
cancermine <- read.table("data/cancermine_genelist_results.tsv", header=T, sep="\t")
cancermine_TSG<-as.character(cancermine[cancermine$tumor_suppressor_citations>0,1])
cancermine_OG<-as.character(cancermine[cancermine$oncogene_citations>0,1])

TSG_sets=list("DORGE-TSGs"=as.character(new_predicted_tsgs[,1]), "CGC-TSGs"=as.character(known_tsg_genes[,1]), "CancerMine TSGs"=cancermine_TSG,"TSGene TSGs"=as.character(unlist(TSGene)))
temp <- venn.diagram(TSG_sets, fill=c("#984EA3","#377EB8","orange","brown"), alpha=c(0.7,0.7,0.7,0.7), cex=	1, cat.fontface=4, filename=NULL)


pdf(file="Raw_figures/Figure_2D_DORGE_TSG_cancermine_TSGene_venn.pdf", family="ArialMT", width=3, height=3)
    grid.draw(temp)
dev.off()


###### Figure 2E: DORGE-OGs, Cancermine and ONGene ######
OG_sets=list("DORGE-OGs"=as.character(new_predicted_ogs[,1]), "CGC-OGs"=as.character(known_og_genes[,1]), "CancerMine OGs"=cancermine_OG, "ONGene OGs"=as.character(unlist(ONGene)))

temp <- venn.diagram(OG_sets, fill=c("#4DAF4A","#E41A1C","orange","brown"), alpha=c(0.7,0.7,0.7,0.7), cex=1, cat.fontface=4, filename=NULL)

pdf(file="Raw_figures/Figure_2E_DORGE_OG_cancermine_ONGene_venn.pdf", family="ArialMT", width=3, height=3)
    grid.draw(temp)
dev.off()