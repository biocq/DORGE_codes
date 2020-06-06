GUST_OG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/Data_file_S2/GUST_OG.txt", header=F, sep="\t",fill=TRUE,quote = "")
GUST_OG<-as.character(GUST_OG$V1)
GUST_OG_nonredundant<-unique(GUST_OG)
write.table(GUST_OG_nonredundant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/nonredundant_GUST_OG.txt",sep="\t",quote=F,row.names=F,col.names=F)

GUST_TSG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/GUST_TSG.txt", header=F, sep="\t",fill=TRUE,quote = "")
GUST_TSG<-as.character(GUST_TSG$V1)
GUST_TSG_nonredundant<-unique(GUST_TSG)
write.table(GUST_TSG_nonredundant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/nonredundant_GUST_TSG.txt",sep="\t",quote=F,row.names=F,col.names=F)

GUST_driver <-c(GUST_TSG_nonredundant,GUST_OG_nonredundant)
GUST_driver_nonredundant <- unique(GUST_driver)
write.table(GUST_driver_nonredundant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Table_1/prediction/nonredundant_GUST_driver.txt",sep="\t",quote=F,row.names=F,col.names=F)


############################

TUSON_TSG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/TUSON_TSG_reported_by_2020plus.txt", header=F, sep="\t",fill=TRUE,quote = "")
TUSON_TSG<-as.character(TUSON_TSG$V1)
TUSON_TSG_nonredundant<-unique(TUSON_TSG)

TUSON_OG <- read.table("/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/TUSON_OG_reported_by_2020plus.txt", header=F, sep="\t",fill=TRUE,quote = "")
TUSON_OG<-as.character(TUSON_OG$V1)
TUSON_OG_nonredundant<-unique(TUSON_OG)

TUSON_driver <-c(TUSON_TSG_nonredundant,TUSON_OG_nonredundant)
TUSON_driver_nonredundant <- unique(TUSON_driver)
write.table(TUSON_driver_nonredundant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Table_1/prediction/TUSON_driver_reported_by_2020plus.txt",sep="\t",quote=F,row.names=F,col.names=F)


############################

DORGE_TSG <- read.table("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/prediction/All_DORGE_predicted_TSGs.txt", header=F, sep="\t",fill=TRUE,quote = "")
DORGE_TSG<-as.character(DORGE_TSG$V1)

DORGE_OG <- read.table("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/prediction/All_DORGE_predicted_OGs.txt", header=F, sep="\t",fill=TRUE,quote = "")
DORGE_OG<-as.character(DORGE_OG$V1)
DORGE_driver_nonredundant <- unique(c(DORGE_TSG,DORGE_OG))
write.table(DORGE_driver_nonredundant,"/Users/jlyu/Box Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Table_1/prediction/DORGE_predicted_drivers.txt",sep="\t",quote=F,row.names=F,col.names=F)
