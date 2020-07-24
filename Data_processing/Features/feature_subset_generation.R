#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Features/")
profile_all<-read.table("All_features.csv",header=T,sep=",");
subset<-c(1,2:23,25:33,35:37,39:40,42,44,48,50,52:53,55:66,68,69,71,72,74,75,77,78,80,81,83,84,86,87,89,90,92,93,95,96,98)
write.table(profile_all[,subset],"All_features_training.csv",sep=",",quote=F,row.names=F,col.names=T)
