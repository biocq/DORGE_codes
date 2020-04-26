################################################
###### Training feature table generation  ######
################################################

###### Input: All_features.csv: All the candidate features including both training features and also those not used due to high feature correlation.
###### Output: All_training_features.csv: The features used in training DORGE.

options(warn=-1)
################################### Data file S2 ###################################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/Figure_codes/Tables/Data_file_S1/");
profile_all<-read.table("../../All_features.csv",header=T,sep=",");
subset<-c(1,2:23,25:33,35:37,39:40,42,44,48,50,52:53,55:66,68,69,71,72,74,75,77,78,80,81,83,84,86,87,89,90,92,93,95,96,98)
write.table(profile_all[,subset],"All_training_features.csv",sep=",",quote=F,row.names=F,col.names=T)