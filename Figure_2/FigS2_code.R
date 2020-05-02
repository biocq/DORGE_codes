#######################################
###### Scatterplots and boxplots ######
#######################################

###### DORGE predictions are compared with TUSON and 20/20+ based on CGC-core and CGC-all genes.


###### Input: ../Gene_set_new.txt: Gene annotation file
###### Input: ../DORGE_prediction.txt: DORGE prediction results
###### Input: data/TUSON_output.txt: TUSON prediction based on the same mutation data in this study
###### Input: data/2020plus_pvals.txt: 20/20+ prediction downloaded from the original paper

###### Output: Figure_S2A_TUSON_DORGE_Rank_scatter_TSG.pdf: Comparison of CGC-TSG gene rankings for DORGE-TSG and TUSON-TSG
###### Output: Figure_S2C_TUSON_DORGE_Rank_scatter_OG.pdf: Comparison of CGC-OG gene rankings for DORGE-OG and TUSON-OG
###### Output: Figure_S2B_plus2020_DORGE_Rank_scatter_TSG.pdf: Comparison of CGC-TSG gene rankings for DORGE-TSG and 20/20+ -TSG
###### Output: Figure_S2D_plus2020_DORGE_Rank_scatter_OG.pdf: Comparison of CGC-OG gene rankings for DORGE-TSG and 20/20+ -TSG
###### Output: Figure_S2E.pdf: Boxplots comparing the ranking of all CGC-TSGs
###### Output: Figure_S2F.pdf Boxplots comparing the ranking of Core CGC-TSGs
###### Output: Figure_S2G.pdf Boxplots comparing the ranking of all CGC-OGs
###### Output: Figure_S2H.pdf Boxplots comparing the ranking of Core CGC-OGs

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ggpubr","scales","dplyr")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library("dplyr"))
suppressMessages(library("ggpubr"))
suppressMessages(library("scales"))

TSG_threshold<-0.62485 #loose FPR=0.01
OG_threshold<-0.7004394 #loose FPR=0.01

#TSG_threshold<-0.8290429 #strict FPR=0.005
#OG_threshold<-0.8679444 #strict FPR=0.005

##################### Figure S2: Scatter plot of ranks from different methods #####################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Figure_2");

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

score <- read.table("../DORGE_prediction.txt", header=T, sep="\t",fill=TRUE,quote = "")
TSG_CGC<-score[,"TSG_all"]
OG_CGC<-score[,"OG_all"]
NG<-score[,"NG"]

TSG_CGC_core<-score[,"TSG_core"]
OG_CGC_core<-score[,"OG_core"]

pvals <- data.frame("Gene"=as.character(score[,"Gene"]),"TSG_probability"=score[,"TSG_probability"])
sorted <- sort(pvals$TSG_probability,decreasing=T,index.return=T)
rank<-1:length(pvals[,1])
Rank_DORGE_TSG<-rep(0,length(pvals[,"Gene"]))
Rank_DORGE_TSG[sorted$ix]<-rank
pvals<-cbind(pvals,Rank_DORGE_TSG)
pvals$Gene<-as.character(pvals$Gene)
pvals2 <- data.frame("Gene"=score[,"Gene"],"OG_probability"=score[,"OG_probability"])
sorted <- sort(pvals2$OG_probability,decreasing=T,index.return=T)
rank<-1:length(pvals2[,1])
Rank_DORGE_OG<-rep(0,length(pvals2[,"Gene"]))
Rank_DORGE_OG[sorted$ix]<-rank
Rank_DORGE_OG<-as.numeric(Rank_DORGE_OG)
pvals2<-cbind(pvals2,Rank_DORGE_OG)
pvals<-cbind(pvals,pvals2[,2:3])

Tuson <- read.table("data/TUSON_output.txt", header=T, sep="\t")
TUSON_final <-data.frame("Gene"=Tuson$Gene,"pval_TSG_TUSON"=Tuson$TUSON.combined.pvalue.TSG,"qval_TSG_TUSON"=Tuson$TUSON.combined.qvalue.TSG,"pval_OG_TUSON"=Tuson$TUSON.combined.pvalue.OG,"qval_OG_TUSON"=Tuson$TUSON.combined.qvalue.OG)
TUSON_final$Gene <- toupper(TUSON_final$Gene)

TUSON_final$pval_TSG_TUSON[is.na(TUSON_final$pval_TSG_TUSON)] <- 1
TUSON_final$pval_OG_TUSON[is.na(TUSON_final$pval_OG_TUSON)] <- 1
join<-left_join(pvals, TUSON_final,by="Gene")
sorted <- sort(join$pval_TSG_TUSON,decreasing=F,index.return=T,na.last=T)
rank<-1:length(pvals[,1])
Rank_TUSON_TSG<-rep(0,length(join[,1]))
Rank_TUSON_TSG[sorted$ix]<-rank
join<-cbind(join,Rank_TUSON_TSG)
sorted <- sort(join$pval_OG_TUSON,decreasing=F,index.return=T,na.last=T)
rank<-1:length(pvals[,1])
Rank_TUSON_OG<-rep(0,length(join[,1]))
Rank_TUSON_OG[sorted$ix]<-rank
join<-cbind(join,Rank_TUSON_OG)

category_TSG<-rep("",nrow(score))

category_OG<-rep("",nrow(score))
category_OG[OG_CGC=="1" & TSG_CGC!="1"]<-"CGC-OG"
category_TSG[TSG_CGC=="1" & OG_CGC!="1"]<-"CGC-TSG"
category_OG[OG_CGC_core=="1" & TSG_CGC_core!="1"]<-"Core CGC-OG"
category_TSG[TSG_CGC_core=="1" & OG_CGC_core!="1"]<-"Core CGC-TSG"

index_pTSG<-which(score$TSG_probability>TSG_threshold & TSG_CGC!="1")
#category_TSG[index_pTSG]<-"Novel DORGE-TSG"
category_TSG[category_TSG==""] <-"Other"
index_pOG<-which(score$OG_probability> OG_threshold & OG_CGC!="1")
#category_OG[index_pOG]<-"Novel DORGE-OG"
category_OG[category_OG==""] <-"Other"
cat<-cbind(category_TSG,category_OG)

join<-cbind(join,cat,by="Gene")

num_TSG_prediction<-length(which(score$TSG_probability>TSG_threshold))
num_OG_prediction<-length(which(score$OG_probability>OG_threshold))

########## Figure_S2A: Scatter plot of TSG ranks by DORGE and TUSON ###########
join2<-join[join$category_TSG!="Other",]
Pvalue_correlation_DORGE_vs_TUSON<-cor.test( ~ Rank_DORGE_TSG + Rank_TUSON_TSG, data=join2,method = "spearman",continuity = FALSE,exact = TRUE,conf.level = 0.95)
Pval_corr<-Pvalue_correlation_DORGE_vs_TUSON$p.value
Rho_corr<-Pvalue_correlation_DORGE_vs_TUSON$estimate

pdf("Raw_figures/Figure_S2A_TUSON_DORGE_Rank_scatter_TSG.pdf", family="ArialMT", width=6.3, height=5)
p<-ggscatter(join2, x = "Rank_DORGE_TSG", y = "Rank_TUSON_TSG",color = "category_TSG",size = 0.6,font.label = c(7,"plain"),xlab="DORGE TSG Rank", ylab="TUSON TSG Rank",palette=c("#377EB8","black","#984EA3"),add = "reg.line",add.params = list(linetype = "dashed"),alpha=0.7)+ scale_x_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_TSG_prediction,2000,5000,10000,20000))+ scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_TSG_prediction,2000,5000,10000,20000))+ geom_hline(yintercept=num_TSG_prediction, linetype="dashed", color = "grey")+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
suppressMessages(print(ggpar(p,legend.title = "",legend=c(0.84,0.08),font.legend = c(9, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=8000, y=20, size=4, label=paste("Rho = ", round(Rho_corr,3), "\nP = ",round(Pval_corr,16), sep=" "), color="black",angle = 0)))
garbage <- dev.off()

########## Figure_S2C: Scatter plot of OG ranks by DORGE and TUSON ###########
join2<-join[join$category_OG!="Other",]
Pvalue_correlation_DORGE_vs_TUSON<-cor.test( ~ Rank_DORGE_OG + Rank_TUSON_OG, data=join2,method = "spearman",continuity = FALSE,exact = TRUE,conf.level = 0.95)
Pval_corr<-Pvalue_correlation_DORGE_vs_TUSON$p.value
Rho_corr<-Pvalue_correlation_DORGE_vs_TUSON$estimate

pdf("Raw_figures/Figure_S2C_TUSON_DORGE_Rank_scatter_OG.pdf", family="ArialMT", width=6.3, height=5)
p<-ggscatter(join2, x = "Rank_DORGE_OG", y = "Rank_TUSON_OG",color = "category_OG",size = 0.6,font.label = c(7,"plain"),xlab="DORGE OG Rank", ylab="TUSON OG Rank",palette=c("#E41A1C","black","#4DAF4A"),add = "reg.line",add.params = list(linetype = "dashed"),alpha=0.7)+ scale_x_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_OG_prediction,2000,5000,10000,20000))+ scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_OG_prediction,2000,5000,10000,20000))+ geom_hline(yintercept=num_OG_prediction, linetype="dashed", color = "grey")+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
suppressMessages(print(ggpar(p,legend.title = "",legend=c(0.84,0.08),font.legend = c(9, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=8000, y=20, size=4, label=paste("Rho = ", round(Rho_corr,3), "\nP = ",round(Pval_corr,7), sep=" "), color="black",angle = 0)))
garbage <- dev.off()

###############################################################################

plus2020 <- read.table("data/2020plus_pvals.txt", header=T, sep="\t")
plus2020_final <-data.frame("Gene"=plus2020$Gene,"pval_TSG_2020plus"=plus2020$tsg_p_value,"qval_TSG_2020plus"=plus2020$tsg_q_value,"pval_OG_2020plus"=plus2020$oncogene_p_value,"qval_OG_2020plus"=plus2020$oncogene_q_value)
plus2020_final$Gene <- toupper(plus2020_final$Gene)

plus2020_final$pval_TSG_2020plus[is.na(plus2020_final$pval_TSG_2020plus)] <- 1
plus2020_final$pval_OG_2020plus[is.na(plus2020_final$pval_OG_2020plus)] <- 1

join<-left_join(join, plus2020_final,by="Gene")

sorted <- sort(join$pval_TSG_2020plus,decreasing=F,index.return=T,na.last=T)
rank<-1:length(pvals[,1])
Rank_2020plus_TSG<-rep(0,length(join[,1]))
Rank_2020plus_TSG[sorted$ix]<-rank
Rank_2020plus_TSG<-as.numeric(Rank_2020plus_TSG)
join<-cbind(join,Rank_2020plus_TSG)

sorted <- sort(join$pval_OG_2020plus,decreasing=F,index.return=T,na.last=T)
rank<-1:length(pvals[,1])
Rank_2020plus_OG<-rep(0,length(join[,1]))
Rank_2020plus_OG[sorted$ix]<-rank
Rank_2020plus_OG<-as.numeric(Rank_2020plus_OG)
join<-cbind(join,Rank_2020plus_OG)

########## Figure_S2B: Scatter plot of TSG ranks by DORGE and 20/20+ ###########
join2<-join[join$category_TSG!="Other",]
Pvalue_correlation_DORGE_vs_2020plus<-cor.test( ~ Rank_DORGE_TSG + Rank_2020plus_TSG, data=join2,method = "spearman",continuity = FALSE,exact = TRUE,conf.level = 0.95)
Pval_corr<-Pvalue_correlation_DORGE_vs_2020plus$p.value
Rho_corr<-Pvalue_correlation_DORGE_vs_2020plus$estimate

pdf("Raw_figures/Figure_S2B_plus2020_DORGE_Rank_scatter_TSG.pdf", family="ArialMT", width=6.3, height=5)
p<-ggscatter(join2, x = "Rank_DORGE_TSG", y = "Rank_2020plus_TSG",color = "category_TSG",size = 0.6,font.label = c(7,"plain"),xlab="DORGE TSG Rank", ylab="20/20+ TSG Rank",palette=c("#377EB8","black","#984EA3"),add = "reg.line",add.params = list(linetype = "dashed"),alpha=0.7)+ scale_x_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_TSG_prediction,2000,5000,10000,20000))+ scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_TSG_prediction,2000,5000,10000,20000))+ geom_hline(yintercept=num_TSG_prediction, linetype="dashed", color = "grey")+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
suppressMessages(print(ggpar(p,legend.title = "",legend=c(0.84,0.08),font.legend = c(9, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=8000, y=20, size=4, label=paste("Rho = ", round(Rho_corr,3), "\nP = ",round(Pval_corr,22), sep=" "), color="black",angle = 0)))
garbage <- dev.off()

########## Figure_S2D: Scatter plot of OG ranks by DORGE and 20/20+ ###########
join2<-join[join$category_OG!="Other",]
Pvalue_correlation_DORGE_vs_2020plus<-cor.test( ~ Rank_DORGE_OG + Rank_2020plus_OG, data=join2,method = "spearman",continuity = FALSE,exact = TRUE,conf.level = 0.95)
Pval_corr<-Pvalue_correlation_DORGE_vs_2020plus$p.value
Rho_corr<-Pvalue_correlation_DORGE_vs_2020plus$estimate

pdf("Raw_figures/Figure_S2D_plus2020_DORGE_Rank_scatter_OG.pdf", family="ArialMT", width=6.3, height=5)
p<-ggscatter(join2, x = "Rank_DORGE_OG", y = "Rank_2020plus_OG",color = "category_OG",size = 0.6,font.label = c(7,"plain"),xlab="DORGE OG Rank", ylab="20/20+ OG Rank",palette=c("#E41A1C","black","#4DAF4A"),add = "reg.line",add.params = list(linetype = "dashed"),alpha=0.7)+ scale_x_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_OG_prediction,2000,5000,10000,20000))+ scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_OG_prediction,2000,5000,10000,20000))+ geom_hline(yintercept=num_OG_prediction, linetype="dashed", color = "grey")+ border()+ theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),legend.background = element_rect(colour = "black",fill = "transparent",size = 0.3, linetype = "solid"),legend.key.size = unit(0.4, "cm"),legend.margin = margin(0.1, 0.1, 0.1, 0.1))
suppressMessages(print(ggpar(p,legend.title = "",legend=c(0.84,0.08),font.legend = c(9, "plain", "black"))+ rremove("legend.title")+ annotate(geom="text", x=8000, y=20, size=4, label=paste("Rho = ", round(Rho_corr,3), "\nP = ",round(Pval_corr,8), sep=" "), color="black",angle = 0)))
garbage <- dev.off()

########## Figure_S2E: Comparison of tools on all CGC-TSGs ###########
inverse = function(yt){
        y <- ((abs(yt) * 0.08 + 1)  ^ (1 / 0.08) - 1) * sign(yt)
        return(median(y))
}

join2<-join[join$category_TSG=="CGC-TSG",]

ranks<-c(join2$Rank_DORGE_TSG,join2$Rank_2020plus_TSG,join2$Rank_TUSON_TSG)
categories<-c(rep("DORGE_TSG",length(join2$Rank_DORGE_TSG)),rep("Plus2020_TSG",length(join2$Rank_2020plus_TSG)),rep("TUSON_TSG",length(join2$Rank_TUSON_TSG)))
combined<-data.frame("Ranks"=ranks,"Cat"=categories)

pdf("Raw_figures/Figure_S2E_comparison_on_all_CGC_TSGs.pdf", family="ArialMT", width=2.5, height=3)
my_comparisons <- list(c("DORGE_TSG", "TUSON_TSG"),c("DORGE_TSG","Plus2020_TSG"))
p<-ggboxplot(data = combined,x = "Cat", y="Ranks",color="black",outlier.size=0.2,width=0.8,fill="Cat",lwd=0.25)+ labs(x = "", y = "Rank")+ scale_x_discrete(labels= c("DORGE TSG","20/20+ TSG","TUSON TSG")) + scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_TSG_prediction,2000,5000,10000,20000))+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 0, hjust = 1, colour = "black"),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "l"),comparisons = my_comparisons,size=2)+stat_summary(fun.data = function(x) data.frame(y=15, size=2, label = paste("Median = ",round(inverse(x),1))), geom="text")+ ggtitle("CGC-TSG")
suppressMessages(print(ggpar(p,legend.title = "")+ rremove("legend")))
garbage <- dev.off()

########## Figure_S2F: Comparison of tools on Core CGC-TSGs ###########
inverse = function(yt){
        y <- ((abs(yt) * 0.08 + 1)  ^ (1 / 0.08) - 1) * sign(yt)
        return(median(y))
}

join2<-join[join$category_TSG=="Core CGC-TSG",]

ranks<-c(join2$Rank_DORGE_TSG,join2$Rank_2020plus_TSG,join2$Rank_TUSON_TSG)
categories<-c(rep("DORGE_TSG",length(join2$Rank_DORGE_TSG)),rep("Plus2020_TSG",length(join2$Rank_2020plus_TSG)),rep("TUSON_TSG",length(join2$Rank_TUSON_TSG)))
combined<-data.frame("Ranks"=ranks,"Cat"=categories)
my_comparisons <- list(c("DORGE_TSG", "TUSON_TSG"),c("DORGE_TSG","Plus2020_TSG"))
pdf("Raw_figures/Figure_S2F_comparison_on_Core_CGC_TSGs.pdf", family="ArialMT", width=2.5, height=3)
p<-ggboxplot(data = combined,x = "Cat", y="Ranks",color="black",outlier.size=0.2,width=0.8,fill="Cat",lwd=0.25)+ labs(x = "", y = "Rank")+ scale_x_discrete(labels= c("DORGE TSG","20/20+ TSG","TUSON TSG")) + scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_TSG_prediction,2000,5000,10000,20000))+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 0, hjust = 1, colour = "black"),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "l"),comparisons = my_comparisons,size=2)+stat_summary(fun.data = function(x) data.frame(y=15, size=2, label = paste("Median=",round(inverse(x),1))), geom="text")+ ggtitle("Core CGC-TSG")
suppressMessages(print(ggpar(p,legend.title = "")+ rremove("legend")))
garbage <- dev.off()

########## Figure_S2G: Comparison of tools on all CGC-OGs ###########
inverse = function(yt){
        y <- ((abs(yt) * 0.08 + 1)  ^ (1 / 0.08) - 1) * sign(yt)
        return(median(y))
}

join2<-join[join$category_OG=="CGC-OG",]

ranks<-c(join2$Rank_DORGE_OG,join2$Rank_2020plus_OG,join2$Rank_TUSON_OG)
categories<-c(rep("DORGE_OG",length(join2$Rank_DORGE_OG)),rep("Plus2020_OG",length(join2$Rank_2020plus_OG)),rep("TUSON_OG",length(join2$Rank_TUSON_OG)))
combined<-data.frame("Ranks"=ranks,"Cat"=categories)
pdf("Raw_figures/Figure_S2G_comparison_on_all_CGC_OGs.pdf", family="ArialMT", width=2.5, height=3)
my_comparisons <- list(c("DORGE_OG", "TUSON_OG"),c("DORGE_OG","Plus2020_OG"))
p<-ggboxplot(data = combined,x = "Cat", y="Ranks",color="black",outlier.size=0.2,width=0.8,fill="Cat",lwd=0.25)+ labs(x = "", y = "Rank")+ scale_x_discrete(labels= c("DORGE OG","20/20+ OG","TUSON OG")) + scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_OG_prediction,2000,5000,10000,20000))+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 0, hjust = 1, colour = "black"),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "l"),comparisons = my_comparisons,size=2)+stat_summary(fun.data = function(x) data.frame(y=15, size=2, label = paste("Median=",round(inverse(x),1))), geom="text")+ ggtitle("CGC-OG")
suppressMessages(print(ggpar(p,legend.title = "")+ rremove("legend")))
garbage <- dev.off()

########## Figure_S2H: Comparison of tools on Core CGC-OGs ###########
inverse = function(yt){
        y <- ((abs(yt) * 0.08 + 1)  ^ (1 / 0.08) - 1) * sign(yt)
        return(median(y))
}

join2<-join[join$category_OG=="Core CGC-OG",]

ranks<-c(join2$Rank_DORGE_OG,join2$Rank_2020plus_OG,join2$Rank_TUSON_OG)
categories<-c(rep("DORGE_OG",length(join2$Rank_DORGE_OG)),rep("Plus2020_OG",length(join2$Rank_2020plus_OG)),rep("TUSON_OG",length(join2$Rank_TUSON_OG)))
combined<-data.frame("Ranks"=ranks,"Cat"=categories)

pdf("Raw_figures/Figure_S2H_comparison_on_Core_CGC_OGs.pdf", family="ArialMT", width=2.5, height=3)
my_comparisons <- list(c("DORGE_OG", "TUSON_OG"),c("DORGE_OG","Plus2020_OG"))
p<-ggboxplot(data = combined,x = "Cat", y="Ranks",color="black",outlier.size=0.2,width=0.8,fill="Cat",lwd=0.25)+ labs(x = "", y = "Rank")+ scale_x_discrete(labels= c("DORGE OG","20/20+ OG","TUSON OG")) + scale_y_continuous(trans = modulus_trans(0.08),breaks = c(0,100,300,num_OG_prediction,2000,5000,10000,20000))+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),axis.text.y = element_text(angle = 0, hjust = 1, colour = "black"),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(fill=FALSE)+stat_compare_means(aes(label = paste0(..p.format..)),method.args = list(alternative = "l"),comparisons = my_comparisons,size=2)+stat_summary(fun.data = function(x) data.frame(y=15, size=2, label = paste("Median=",round(inverse(x),1))), geom="text")+ ggtitle("Core CGC-OG")
suppressMessages(print(ggpar(p,legend.title = "")+ rremove("legend")))
garbage <- dev.off()

#write.table(join,"ranking_comparison.txt",sep="\t",quote=F,row.names=F)