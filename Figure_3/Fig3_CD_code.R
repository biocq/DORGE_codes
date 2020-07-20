######################################
##### Gene enrichment evaluation #####
######################################

###### The gene enrichment analyses for different gene categories to independently evaluate DORGE prediction.

###### Input: ../All_features.csv: Feature table compiled in this study
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/ER_genes.txt: Epigenetic regulator (ER) gene set (Compiled from https://genome.cshlp.org/content/29/4/532.long and https://academic.oup.com/database/article/doi/10.1093/database/bav067/2433200; Available at Evaluation_processing folder DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)
###### Input: data/driver-table.txt: Sleeping Beauty (SB) gene set from Sleeping Beauty Cancer Driver Database (SBCDDB) (https://pubmed.ncbi.nlm.nih.gov/23770567/; Available at Evaluation_processing folder DORGE_pipeline project https://github.com/biocq/DORGE_pipeline)

###### Output: Figure_3C_ER_enrichment.pdf: Enrichment of ERs in different gene categories 
###### Output: Figure_3D_SB_enrichment.pdf: Enrichment of SBs in different gene categories

options(warn=-1)

### Install missing packages

pkgs <-  c("ComplexHeatmap")
installed_pkgs <- installed.packages()

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
}

pkgs <-  c("ggpubr","dplyr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("dplyr"))
suppressMessages(library("ggpubr"))


TSG_threshold<-0.6233374 #FPR=0.01
OG_threshold<-0.6761319 #FPR=0.01

##################### Figure 3C: Epigenetic regulator (ER) gene enrichment #####################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_3");

anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","TSG_all","OG_all","NG")]
index_NG<-which(allgene$NG=="1")
index_pOG<-which(allgene$OG_probability> OG_threshold & allgene$TSG_probability < TSG_threshold & allgene$OG_all!="1")
index_pTSG<-which(allgene$TSG_probability>TSG_threshold & allgene$OG_probability< OG_threshold & allgene$TSG_all!="1")
index_TSG<-which(allgene$TSG_all=="1")
index_OG<-which(allgene$OG_all=="1")
index_driver<-which(allgene$TSG_all=="1" | allgene$OG_all=="1")

CGC_TSG<-as.character(allgene$Gene[index_TSG])
CGC_OG<-as.character(allgene$Gene[index_OG])
DORGE_novel_TSG<-as.character(allgene$Gene[index_pTSG])
DORGE_novel_OG<-as.character(allgene$Gene[index_pOG])
NG<-as.character(allgene$Gene[anno[,"NG"]==1])
driver<-as.character(allgene$Gene[index_driver])

ER_genes <- read.table("data/ER_genes.txt", header=F, sep="\t")
ER_genes <- as.character(ER_genes$V1)

non_CGC_TSG_genes<-setdiff(allgene[,1],CGC_TSG)
non_CGC_OG_genes<-setdiff(allgene[,1],CGC_OG)
non_pTSG_genes<-setdiff(allgene[,1],DORGE_novel_TSG)
non_pOG_genes<-setdiff(allgene[,1],DORGE_novel_OG)
non_driver_genes<-setdiff(allgene[,1],driver)
non_NG_genes<-setdiff(allgene[,1],NG)
non_ER_genes<-setdiff(allgene[,1],ER_genes)

knownTSG_ER<-length(which(CGC_TSG%in%ER_genes))
knownTSG_nonER<-length(which(CGC_TSG%in%non_ER_genes))
nonknownTSG_ER<-length(which(non_CGC_TSG_genes%in%ER_genes))
nonknownTSG_nonER<-length(which(non_CGC_TSG_genes%in%non_ER_genes))
fisher<-matrix(c(knownTSG_ER,knownTSG_nonER,nonknownTSG_ER,nonknownTSG_nonER),nrow = 2,dimnames =list(c("TSG", "nonTSG"),c("ER", "nonER")))
pval_knownTSG<-fisher.test(round(fisher*200/length(index_TSG)), alternative = "g")$p.value

knownOG_ER<-length(which(CGC_OG%in%ER_genes))
knownOG_nonER<-length(which(CGC_OG%in%non_ER_genes))
nonknownOG_ER<-length(which(non_CGC_OG_genes%in%ER_genes))
nonknownOG_nonER<-length(which(non_CGC_OG_genes%in%non_ER_genes))
fisher<-matrix(c(knownOG_ER,knownOG_nonER,nonknownOG_ER,nonknownOG_nonER),nrow = 2,dimnames =list(c("OG", "nonOG"),c("ER", "nonER")))
pval_knownOG<-fisher.test(round(fisher*200/length(index_OG)), alternative = "g")$p.value

novelTSG_ER<-length(which(DORGE_novel_TSG%in%ER_genes))
novelTSG_nonER<-length(which(DORGE_novel_TSG%in%non_ER_genes))
nonnovelTSG_ER<-length(which(non_pTSG_genes%in%ER_genes))
nonnovelTSG_nonER<-length(which(non_pTSG_genes%in%non_ER_genes))
fisher<-matrix(c(novelTSG_ER,novelTSG_nonER,nonnovelTSG_ER,nonnovelTSG_nonER),nrow = 2,dimnames =list(c("pTSG", "nonpTSG"),c("ER", "nonER")))
pval_novelTSG<-fisher.test(round(fisher*200/length(index_pTSG)), alternative = "g")$p.value

novelOG_ER<-length(which(DORGE_novel_OG%in%ER_genes))
novelOG_nonER<-length(which(DORGE_novel_OG%in%non_ER_genes))
nonnovelOG_ER<-length(which(non_pOG_genes%in%ER_genes))
nonnovelOG_nonER<-length(which(non_pOG_genes%in%non_ER_genes))
fisher<-matrix(c(novelOG_ER,novelOG_nonER,nonnovelOG_ER,nonnovelOG_nonER),nrow = 2,dimnames =list(c("pOG", "nonpOG"),c("ER", "nonER")))
pval_novelOG<-fisher.test(round(fisher*200/length(index_pOG)), alternative = "g")$p.value

driver_ER<-length(which(driver%in%ER_genes))
driver_nonER<-length(which(driver%in%non_ER_genes))
nondriver_ER<-length(which(non_driver_genes%in%ER_genes))
nondriver_nonER<-length(which(non_driver_genes%in%non_ER_genes))
fisher<-matrix(c(driver_ER,driver_nonER,nondriver_ER,nondriver_nonER),nrow = 2,dimnames =list(c("driver", "nondriver"),c("ER", "nonER")))
pval_driver<-fisher.test(round(fisher*200/length(index_driver)), alternative = "g")$p.value

NG_ER<-length(which(driver%in%NG))
NG_nonER<-length(which(driver%in%non_NG_genes))
nonNG_ER<-length(which(non_driver_genes%in%NG))
nonNG_nonER<-length(which(non_driver_genes%in%non_NG_genes))
fisher<-matrix(c(NG_ER,NG_nonER,nonNG_ER,nonNG_nonER),nrow = 2,dimnames =list(c("NG", "nonNG"),c("ER", "nonER")))
pval_NG<-fisher.test(round(fisher*200/length(index_NG)), alternative = "g")$p.value

pval<-c(pval_novelTSG,pval_novelOG,pval_NG/2,pval_driver,pval_knownTSG,pval_knownOG)
pval<- -log10(pval)
id<-c("Novel DORGE-TSG","Novel DORGE-OG","NG","CGC driver","CGC-TSG","CGC-OG")
dat<-data.frame(id,pval)
pdf("Raw_figures/Figure_3C_ER_enrichment.pdf", family="ArialMT", width=4, height=3)
ggplot(data=dat, aes(x=as.factor(id), y=pval,fill=as.factor(id)))+geom_bar(stat="identity", alpha=1, position="dodge") + coord_flip()+labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y=element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ scale_x_discrete(limits = rev(id))+scale_fill_manual(values=c("#333333","#E41A1C","#377EB8","#CCCCCC","#4DAF4A","#984EA3"))+ ggtitle("ER enrichment")+ annotate("text", x = 1, y = 10,label = paste("n = ",knownOG_ER,sep=""), parse = F,size =3)+annotate("text", x = 2, y = 16,color="white",label = paste("n = ",knownTSG_ER,sep=""), parse = F,size =3)+annotate("text", x = 3, y = 9,color="white",label = paste("n = ",driver_ER,sep=""), parse = F,size =3)+annotate("text", x = 4, y = 3,label = paste("n = ",NG_ER,sep=""), parse = F,size =3)+annotate("text", x = 5, y = 5,label = paste("n = ",novelOG_ER,sep=""), parse = F,size =3)+annotate("text", x = 6, y = 9,label = paste("n = ",novelTSG_ER,sep=""), parse = F,size =3)
garbage <- dev.off()


##################### Figure 3D: Sleeping Beauty (SB) gene enrichment #####################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_3");

SBs <- read.table("data/driver-table.txt", header=T, sep="\t",fill=TRUE,quote = "")
SBs <- mutate_each(SBs, list(toupper))
SBs2 <- SBs[SBs$pattern=="INACTIVATING",]
inactivated_genes<-as.character(SBs2$gene)

anno <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")

allgene<-anno[,c("Gene","TSG_probability","OG_probability","TSG_core","OG_core","TSG_all","OG_all","NG")]

index_TSG<-which(allgene$TSG_all=="1")# & allgene$OG_all!="1"
index_NG<-which(allgene$NG=="1")
index_pTSG<-which(allgene$TSG_probability>TSG_threshold & allgene$OG_probability< OG_threshold & allgene$TSG_all!="1")

CGC_TSG<-as.character(anno[index_TSG,1])
NGs<-as.character(anno[index_NG,1])
DORGE_novel_TSGs<-as.character(anno[index_pTSG,1])

non_CGC_TSG_genes<-setdiff(allgene[,1],CGC_TSG)
non_pTSG_genes<-setdiff(allgene[,1],DORGE_novel_TSGs)
non_NG_genes<-setdiff(allgene[,1],NGs)
noninactivated_genes<-setdiff(allgene[,1],inactivated_genes)

knownTSG_SB<-length(which(CGC_TSG%in%inactivated_genes))
knownTSG_nonSB<-length(which(CGC_TSG%in%noninactivated_genes))
nonknownTSG_SB<-length(which(non_CGC_TSG_genes%in%inactivated_genes))
nonknownTSG_nonSB<-length(which(non_CGC_TSG_genes%in%noninactivated_genes))
fisher<-matrix(c(knownTSG_SB,knownTSG_nonSB,nonknownTSG_SB,nonknownTSG_nonSB),nrow = 2,dimnames =list(c("TSG", "nonTSG"),c("SB", "nonSB")))
pval_knownTSG<-fisher.test(round(fisher*200/length(index_TSG)), alternative = "g")$p.value

novelTSG_SB<-length(which(DORGE_novel_TSGs%in%inactivated_genes))
novelTSG_nonSB<-length(which(DORGE_novel_TSGs%in%noninactivated_genes))
nonnovelTSG_SB<-length(which(non_pTSG_genes%in%inactivated_genes))
nonnovelTSG_nonSB<-length(which(non_pTSG_genes%in%noninactivated_genes))
fisher<-matrix(c(novelTSG_SB,novelTSG_nonSB,nonnovelTSG_SB,nonnovelTSG_nonSB),nrow = 2,dimnames =list(c("TSG", "nonTSG"),c("SB", "nonSB")))
pval_novelTSG<-fisher.test(round(fisher*200/length(index_pTSG)), alternative = "g")$p.value

NG_SB<-length(which(NGs%in%inactivated_genes))
NG_nonSB<-length(which(NGs%in%noninactivated_genes))
nonNG_SB<-length(which(non_NG_genes%in%inactivated_genes))
nonNG_nonSB<-length(which(non_NG_genes%in%noninactivated_genes))
fisher<-matrix(c(NG_SB,NG_nonSB,nonNG_SB,nonNG_nonSB),nrow = 2,dimnames =list(c("NG", "nonNG"),c("SB", "nonSB")))
pval_NG<-fisher.test(round(fisher*200/length(index_NG)), alternative = "g")$p.value/2

pval<-c(pval_knownTSG,pval_NG,pval_novelTSG)
pval<- -log10(pval)
id<-c("CGC-TSG","NG","Novel DORGE-TSG")
dat<-data.frame(id,pval)
dat$id<-factor(dat$id,levels=c("CGC-TSG","Novel DORGE-TSG","NG"))

pdf("Raw_figures/Figure_3D_SB_enrichment.pdf",family="ArialMT", width=1.5, height=3)
ggplot(data=dat, aes(x=as.factor(id), y=pval,fill=as.factor(id)))+geom_bar(stat="identity", alpha=1, position="dodge") + labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ theme_bw() + theme(plot.title = element_text(size=10),axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"),legend.position = "none",panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ scale_x_discrete(limits = id)+scale_fill_manual(values=c("#377EB8","#984EA3","#CCCCCC"))+ ggtitle("SB enrichment")+ annotate("text", x = 1,y = 19.5,color="black",label = paste("n = ",knownTSG_SB,sep=""), parse = F,size =2)+annotate("text", x = 2, y = 1.5,color="black",label = paste("n = ",NG_SB,sep=""), parse = F,size =2)+annotate("text", x = 3, y = 24.5,color="black",label = paste("n = ",novelTSG_SB,sep=""),size=2, parse = F)
garbage <- dev.off()