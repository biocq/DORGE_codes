###############################
######### Boxplots ############
###############################

###### The barplots of top feature groups for TSG and OG predictions are presented.
###### Input: TSG_top4_feature_groups.txt
###### Input: OG_top7_feature_groups.txt
###### Output: Figure_1A_TSG_feature_groups.pdf
###### Output: Figure_1B_OG_feature_groups.pdf

options(warn=-1)

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("plyr","ggpubr","cowplot")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs), repos = "http://cran.us.r-project.org")
}

suppressMessages(library(ggpubr));
suppressMessages(library(plyr));
suppressMessages(library(cowplot));

################################## Figure 1A Feature groups contributing most to the TSG prediction ###############################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Figure_1");

TSG_features<-read.table("data/TSG_top3_feature_groups.txt",header=T,sep="\t");
TSG_features[,1]<-gsub("_"," ",TSG_features[,1])
TSG_features$Group <- as.factor(TSG_features$Group)
TSG_features$Group<-factor(TSG_features$Group,levels=c("Group 1", "Group 2", "Group 3"))
TSG_features$minuslogP<-rep(0,nrow(TSG_features))
TSG_features$minuslogP<- -log10(TSG_features$Wilcoxon_Test_P_value)
TSG_features$hm<-rep("",nrow(TSG_features))
TSG_features$hm[grep("length", TSG_features$Feature)]<-"length"
TSG_features$hm[grep("Height", TSG_features$Feature)]<-"height"
TSG_features$short_name<-TSG_features$Feature
TSG_features$short_name<-gsub(" peaks", "", TSG_features$short_name)
TSG_features$short_name<-gsub("Height of ", "", TSG_features$short_name)
TSG_features$short_name<-gsub(" peak length", "", TSG_features$short_name)
TSG_features$short_name<-factor(TSG_features$short_name,levels=(c("H3K79me2","H3K36me3","H4K20me1","H3K4me2","H3K4me1","H3K9ac","H3K27ac","H3K4me3")))
#group_values<-data.frame(group=c("Group 1","Group 2","Group 3","Group 4"),values=c(0.07106825,0.008291412,0.00680956,0.005146332))
#group_values$group<-factor(group_values$group,levels=rev(c("Group 1","Group 2","Group 3","Group 4")))

p01<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 2, label = "Group 1:\nAUPRC\nreduction = 0.0721"),aes(x, y, label = label),hjust = 0.0, vjust = 0.5, size = 23,color = "black",inherit.aes = FALSE)
p02<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 2, label = "Group 2:\nAUPRC\nreduction = 0.0080"),aes(x, y, label = label),hjust = 0.0, vjust = 0.5, size = 23,color = "black",inherit.aes = FALSE)
p04<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 2, label = "Group 3:\nAUPRC\nreduction = 0.0054"),aes(x, y, label = label),hjust = 0.0, vjust = 0.5, size = 23,color = "black",inherit.aes = FALSE)

p1<-ggbarplot(TSG_features[TSG_features$hm=="length",],x="short_name",y="minuslogP",sort.val = "desc",fill = "darkgreen",color="transparent",sort.by.groups = TRUE,width = 0.7)+ theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="Peak length",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ rremove("x.ticks")+ rremove("x.text")+ggtitle("Histone modification")

p2<-ggbarplot(TSG_features[TSG_features$hm=="height",],x="short_name",y="minuslogP",sort.val = "none",fill = "blue",color="transparent",sort.by.groups = F,width = 0.7)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="Peak height",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")

p4<-ggbarplot(TSG_features[TSG_features$Group=="Group 2",],x="Feature",y="minuslogP",sort.val = "desc",fill = "purple",color="transparent",sort.by.groups = TRUE,width = 0.6)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ggtitle("Phenotype")


p5<-ggbarplot(TSG_features[TSG_features$Group=="Group 3",],x="Feature",y="minuslogP",sort.val = "desc",fill = "darkred",color="transparent",sort.by.groups = TRUE,width = 0.7)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ggtitle("Missense mutations")

combined<-plot_grid(
	NULL,
	plot_grid(
		plot_grid(
			plot_grid(p1, p2, labels = "", align = 'v', rel_heights=c(1,1.0),nrow=2)
			,NULL,ncol=2, align = 'h'
		), 
		p01,
		plot_grid(p4, NULL,ncol=2, align = 'v',rel_widths=c(1,1.3)),
		p02,
		plot_grid(p5, NULL,ncol=2, align = 'v',rel_widths=c(1,1.3)),
		p04,
		axis="l", labels = "", align = 'v', rel_heights=c(0.0065,0.008,0.006,0.005,0.006,0.005),nrow=6
	),
	NULL, labels = "", rel_widths = c(3.5,4,4),ncol=3
)

save_plot("Raw_figures/Figure_1A_TSG_feature_groups.pdf", combined, family="ArialMT", ncol = 2, nrow = 3, base_height = 20,limitsize = FALSE)

################################## Figure 1B Feature groups contributing most to the OG prediction ###############################


OG_features<-read.table("data/OG_top5_feature_groups.txt",header=T,sep="\t");
OG_features[,1]<-gsub("_"," ",OG_features[,1])
OG_features$Group <- as.factor(OG_features$Group)
OG_features$Group<-factor(OG_features$Group,levels=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"))
OG_features$minuslogP<-rep(0,nrow(OG_features))
OG_features$minuslogP<- -log10(OG_features$Wilcoxon_Test_P_value)
OG_features$hm<-rep("",nrow(OG_features))
OG_features$hm[grep("length", OG_features$Feature)]<-"length"
OG_features$hm[OG_features$Feature=="log gene length"]<-""
OG_features$hm[OG_features$Feature=="log CDS length"]<-""
OG_features$hm[grep("Height", OG_features$Feature)]<-"height"
OG_features$short_name<-OG_features$Feature
OG_features$short_name<-gsub(" peaks", "", OG_features$short_name)
OG_features$short_name<-gsub("Height of ", "", OG_features$short_name)
OG_features$short_name<-gsub(" peak length", "", OG_features$short_name)
OG_features$short_name<-factor(OG_features$short_name,levels=(c("H3K79me2","H3K36me3","H4K20me1","H3K27ac","H3K4me1","H3K9ac","H3K4me2","H3K4me3")))

p01<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 1, label = "Group 1:\nAUPRC\nreduction = 0.0509"),aes(x, y, label = label),hjust = 0.0, vjust = 0.5, size = 23,color = "black",inherit.aes = FALSE)
p02<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 1, label = "Group 2:\nAUPRC\nreduction = 0.0316"),aes(x, y, label = label),hjust = 0.0, vjust = 0.0, size = 23,color = "black",inherit.aes = FALSE)
p03<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 1, label = "Group 3:\nAUPRC\nreduction = 0.0231"),aes(x, y, label = label),hjust = 0.0, vjust = 0.0, size = 23,color = "black",inherit.aes = FALSE)
p04<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 1, label = "Group 4:\nAUPRC\nreduction = 0.0123"),aes(x, y, label = label),hjust = 0.0, vjust = 0.0, size = 23,color = "black",inherit.aes = FALSE)
p05<-ggplot() + theme_void()+ geom_text(data = data.frame(x = 1, y = 1, label = "Group 5:\nAUPRC\nreduction = 0.0072"),aes(x, y, label = label),hjust = 0.0, vjust = 0.0, size = 23,color = "black",inherit.aes = FALSE)


p1<-ggbarplot(OG_features[OG_features$Group=="Group 1",],x="Feature",y="minuslogP",sort.val = "desc",fill = "darkred",color="white",sort.by.groups = TRUE,width = 0.8)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.length=unit(0.5, "cm"),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ggtitle("Missense mutations")

p2<-ggbarplot(OG_features[OG_features$Group=="Group 2",],x="Feature",y="minuslogP",sort.val = "desc",fill = "cyan4",color="transparent",sort.by.groups = TRUE,width = 0.75)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.length=unit(0.5, "cm"),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ggtitle("Genomics")

p3<-ggbarplot(OG_features[OG_features$Group=="Group 3",],x="Feature",y="minuslogP",sort.val = "desc",fill = "coral",color="transparent",sort.by.groups = TRUE,width = 0.7)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.length=unit(0.5, "cm"),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ggtitle("Super enhancer")

p4<-ggbarplot(OG_features[OG_features$Group=="Group 4",],x="Feature",y="minuslogP",sort.val = "desc",fill = "orange",color="transparent",sort.by.groups = TRUE,width = 0.5)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.length=unit(0.5, "cm"),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ggtitle("DNA methylation")


p6<-ggbarplot(OG_features[OG_features$hm=="length",],x="short_name",y="minuslogP",sort.val = "desc",fill = "darkgreen",color="transparent",sort.by.groups = TRUE,width = 0.75)+ theme_bw() + theme(axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35), panel.border = element_blank(),panel.grid.major = element_blank(),axis.ticks.length=unit(0.5, "cm"),panel.grid.minor = element_blank(),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="Peak length",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")+ rremove("x.ticks")+ rremove("x.text")+ggtitle("Histone modification")

p7<-ggbarplot(OG_features[OG_features$hm=="height",],x="short_name",y="minuslogP",sort.val = "none",fill = "blue",color="transparent",sort.by.groups = F,width = 0.75)+ theme_bw() + theme(axis.ticks.x = element_line(size = 0.5),axis.ticks.y = element_line(size = 0.5),plot.margin = margin(0.01, 0.01, 0.01, 0.01),text = element_text(size=35),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.length=unit(0.5, "cm"),axis.line = element_line(colour ="black"),axis.text.y = element_text(colour = "black", size = 35),axis.text.x = element_text(angle = 30, size = 35, hjust = 1, colour = "black"))+labs(x="Peak height",y=expression("-log"[10]~italic(P)~"-value"))+ guides(fill=guide_legend(title=""))+ rremove("legend")


combined<-plot_grid(
	NULL,
	plot_grid(
		plot_grid(p1, NULL,ncol=2, align = 'v',rel_widths=c(1,1.4)),
		p01,
		plot_grid(p2, NULL,ncol=2, align = 'v',rel_widths=c(1,0.9)),
		p02,
		plot_grid(p3, NULL,ncol=2, align = 'v',rel_widths=c(1,6.0)),
		p03,
		plot_grid(p4, NULL,ncol=2, align = 'v',rel_widths=c(1,4.5)),
		p04,
		plot_grid(
			plot_grid(p6, p7, labels = "", align = 'hv', rel_heights=c(1,1),nrow=2)
			,NULL,ncol=2, align = 'h',rel_widths=c(1,1.3)
		), 
		p05,
		axis="l", labels = "", align = 'v', rel_heights=c(0.0036,0.005, 0.0036,0.005, 0.0030,0.005, 0.0025,0.004, 0.0036, 0.005),nrow=10
	),
	NULL, labels = "", rel_widths = c(4,4,4.),ncol=3
)


save_plot("Raw_figures/Figure_1B_OG_feature_groups.pdf", combined, family="ArialMT", ncol = 3, nrow = 7, base_height = 25,limitsize = FALSE)