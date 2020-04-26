GUST_OG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/GUST_OG.txt", header=F, sep="\t",fill=TRUE,quote = "")
GUST_OG<-as.character(GUST_OG$V1)
GUST_OG_nonredudant<-unique(GUST_OG)
write.table(GUST_OG_nonredudant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/nonredudant_GUST_OG.txt",sep="\t",quote=F,row.names=F,col.names=F)

GUST_TSG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/GUST_TSG.txt", header=F, sep="\t",fill=TRUE,quote = "")
GUST_TSG<-as.character(GUST_TSG$V1)
GUST_TSG_nonredudant<-unique(GUST_TSG)
write.table(GUST_TSG_nonredudant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/nonredudant_GUST_TSG.txt",sep="\t",quote=F,row.names=F,col.names=F)

GUST_driver <-c(GUST_TSG_nonredudant,GUST_OG_nonredudant)
GUST_driver_nonredudant <- unique(GUST_driver)
write.table(GUST_driver_nonredudant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/nonredudant_GUST_driver.txt",sep="\t",quote=F,row.names=F,col.names=F)


############################

TUSON_TSG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/TUSON_TSG_reported_by_2020plus.txt", header=F, sep="\t",fill=TRUE,quote = "")
TUSON_TSG<-as.character(TUSON_TSG$V1)
TUSON_TSG_nonredudant<-unique(TUSON_TSG)

TUSON_OG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/TUSON_OG_reported_by_2020plus.txt", header=F, sep="\t",fill=TRUE,quote = "")
TUSON_OG<-as.character(TUSON_OG$V1)
TUSON_OG_nonredudant<-unique(TUSON_OG)

TUSON_driver <-c(TUSON_TSG_nonredudant,TUSON_OG_nonredudant)
TUSON_driver_nonredudant <- unique(TUSON_driver)
write.table(TUSON_driver_nonredudant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/cancer_driver_prediction/q0point1/TUSON_driver_reported_by_2020plus.txt1",sep="\t",quote=F,row.names=F,col.names=F)