######################
##### PRC curves #####
######################


###### DORGE predicted genes are characterized by different features.

###### Input: data/TSG.rdata: PRC data for TSG prediction
###### Input: data/OG.rdata: PRC data for OG prediction
###### Output: Figure_2B_PRC_TSGs.pdf: PRC curves for TSGs using all features as well as predefined feature subsets.
###### Output: Figure_2C_PRC_OGs.pdf: PRC curves for OGs using all features as well as predefined feature subsets.

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("PRROC","ggplot2")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

library(PRROC)
library(ggplot2)
################################## Figure 2B: TSG PRC curves ###############################

#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_2");
load("data/TSG.rdata")
prc <- pr.curve(scores.class0=TSG.predictions[[1]], weights.class0=y_TSG, curve=T)
dat <- as.data.frame(prc$curve)
dat<-dat[order(-dat$V1, dat$V2), ]
plot_dat<-data.frame(Recall=dat$V1,Precision=dat$V2,Group=names(TSG.predictions)[1])
	
for (i in 2:length(TSG.predictions)) {
  prc <- pr.curve(scores.class0=TSG.predictions[[i]], weights.class0=y_TSG, curve=T)
  dat <- as.data.frame(prc$curve)
  dat<-dat[ order(-dat$V1, dat$V2), ]
	plot_dat<-rbind(plot_dat,data.frame(Recall=dat$V1,Precision=dat$V2,Group=names(TSG.predictions)[i]))
	#plot(prc, main=names(TSG.predictions)[i])
}

#> names(TSG.predictions)
#"Mutation","Genomics","Phenotype","Epigenetics","wo_Mutation","wo_Genomics","wo_Phenotype","wo_Epigenetics"
#"CRISPR-only","TUSON-TSG","All"

pdf(file="Raw_figures/Figure_2B_PRC_TSGs.pdf", family="ArialMT", width=6.623, height=5)
ggplot(data = plot_dat, aes(x = Recall, y = Precision))+scale_x_continuous(limits = c(0, 1),labels = c("0","0.25", "0.50", "0.75","1"))+scale_y_continuous(limits = c(0, 1),labels = c("0","0.25", "0.50", "0.75","1"))+geom_path(aes(color = Group, linetype = Group))+scale_color_manual(values = c("orange3","red","seagreen3","dodgerblue4","orange3","red","seagreen3","dodgerblue4","darkviolet","gold2","black"))+scale_linetype_manual(values=c("solid","solid","solid","solid","11","11","11","11","1131","solid","solid")) + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_line(size = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)
garbage <- dev.off()

################################## Figure 2C: OG PRC curves ###############################
load("data/OG.rdata")
prc <- pr.curve(scores.class0=OG.predictions[[1]], weights.class0=y_OG, curve=T)
dat <- as.data.frame(prc$curve)
dat<-dat[order(-dat$V1, dat$V2), ]
plot_dat<-data.frame(Recall=dat$V1,Precision=dat$V2,Group=names(OG.predictions)[1])
	
for (i in 2:length(OG.predictions)) {
  prc <- pr.curve(scores.class0=OG.predictions[[i]], weights.class0=y_OG, curve=T)
  dat <- as.data.frame(prc$curve)
  dat<-dat[ order(-dat$V1, dat$V2), ]
	plot_dat<-rbind(plot_dat,data.frame(Recall=dat$V1,Precision=dat$V2,Group=names(OG.predictions)[i]))
	#plot(prc, main=names(OG.predictions)[i])
}

#> names(OG.predictions)
#"Mutation","Genomics","Phenotype","Epigenetics","wo_Mutation","wo_Genomics","wo_Phenotype","wo_Epigenetics"
#"CRISPR-only","TUSON-OG","All"

pdf(file="Raw_figures/Figure_2C_PRC_OGs.pdf", family="ArialMT", width=6.623, height=5)
ggplot(data = plot_dat, aes(x = Recall, y = Precision))+scale_x_continuous(limits = c(0, 1),labels = c("0","0.25", "0.50", "0.75","1"))+scale_y_continuous(limits = c(0, 1),labels = c("0","0.25", "0.50", "0.75","1"))+geom_path(aes(color = Group, linetype = Group))+scale_color_manual(values = c("orange3","red","seagreen3","dodgerblue4","orange3","red","seagreen3","dodgerblue4","darkviolet","gold2","black"))+scale_linetype_manual(values=c("solid","solid","solid","solid","11","11","11","11","1131","solid","solid")) + theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),axis.ticks.x=element_line(size = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)
garbage <- dev.off()