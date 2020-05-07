##########################################################
###### Driver gene prediction evaluation - Table 1  ######
##########################################################

###### Input: gene prediction from prediction folder: Gene lists predicted by different methods in the prediction folder.
###### Output: Table1_driver_gene_performance_on_CGC.txt: Table 1.

options(warn=-1)
################################### Table 1 ###################################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Table_1/");

anno<-read.table("../../Gene_set_new.txt",header=T,sep="\t");
genes<-as.character(anno[,1])
drivers<-genes[which(as.character(anno[,5])==1 | as.character(anno[,6])==1)]

NGs<-genes[which(as.character(anno[,4])==1)]

result<-c()
methods<-c("2020_driver.txt","DORGE_predicted_drivers.txt","ActiveDriver.txt","MutPanning.txt","MuSIC.txt"
,"MutSigCV.txt","OncodriveCLUST.txt","OncodriveFM.txt","OncodriveFML.txt","nonredudant_GUST_driver.txt","TUSON_driver_reported_by_2020plus.txt")


for(i in methods){

		prediction<-read.table(paste("prediction/",i,sep=""),header=F);
		prediction<-as.character(prediction[,1])
		nonprediction<-setdiff(genes,prediction)
		N<-length(prediction)
		correctedly_predicted_driver<-as.character(intersect(prediction,drivers))
		falsely_predicted_driver<-intersect(drivers,nonprediction)
		correctedly_predicted_NG<-intersect(NGs,nonprediction)
		falsely_predicted_NG<-as.character(intersect(NGs,prediction))
				
		#							 Reference	
		#Predicted		driver	  NG
		#driver				A					C
		#nondriver		B					D
		
		A<-length(correctedly_predicted_driver)
		B<-length(falsely_predicted_driver)
		C<-length(falsely_predicted_NG)
		D<-length(correctedly_predicted_NG)
		
		Sensitivity<-round(A/(A+B),3)
		Specificity<-round(D/(C+D),3)
		Precision<-round(A/(A+C),3)
		Accuracy<-round((A+D)/(A+B+C+D),3)
		result<-rbind(result,c(N,Sensitivity,Specificity,Precision,Accuracy))
}
result<-cbind(c("20/20+","DORGE","ActiveDriver","MutPanning","MuSIC","MutSigCV","OncodriveCLUST","OncodriveFM","OncodriveFML","GUST","TUSON"),result)
result<-as.data.frame(result)
colnames(result)<-c("Method","#","Sn","Sp","Precision","Accuracy")
result$Accuracy<-as.numeric(result$Accuracy)
result<-result[order(-result$Accuracy),]

write.table(result,"Table1_driver_gene_performance_on_CGC.txt",sep="\t",quote=F,row.names=F,col.names=T)