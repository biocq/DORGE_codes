###############################################################
###### TSG/OG gene prediction evaluation - Data file S2  ######
###############################################################

###### Input: TSG/OG prediction from prediction folder: Gene lists predicted by different methods in the prediction folder.
###### Input: TSG/OG prediction from topgenes folder: Gene lists predicted by different methods in the topgenes folder, the gene lists are redefined to make the predicted gene numbers the same to DORGE.
###### Output: TSG_gene_performance.txt: TSG prediction performance by different methods on CGC genes.
###### Output: OG_gene_performance.txt: OG prediction performance by different methods on CGC genes.

options(warn=-1)
################################### Data file S2 ###################################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S2/");


anno<-read.table("../../Gene_set_new.txt",header=T,sep="\t");
genes<-as.character(anno[,1])
NGs<-genes[which(as.character(anno[,4])=="1")]
drivers<-genes[which(as.character(anno[,5])=="1")]

result<-c()
methods<-c("All_DORGE_predicted_TSGs.txt","2020_TSG.txt","nonredudant_GUST_TSG.txt","TUSON_TSG_reported_by_2020plus.txt")

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

methods<-c("2020_TSG.txt","TUSON_TSG.txt")

for(i in methods){

		prediction<-read.table(paste("topgenes/",i,sep=""),header=F);
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

result<-cbind(c("DORGE","20/20+","GUST","TUSON","20/20+*","TUSON*"),result)

result<-as.data.frame(result)
colnames(result)<-c("Method","#","Sn","Sp","Precision","Accuracy")
result$Accuracy<-as.numeric(result$Accuracy)
result<-result[c(order(-result[1:4,]$Accuracy),c(6,5)),]

write.table(result,"TSG_gene_performance.txt",sep="\t",quote=F,row.names=F,col.names=T)



drivers<-genes[which(as.character(anno[,6])=="1")]
result<-c()
methods<-c("All_DORGE_predicted_OGs.txt","2020_OG.txt","nonredudant_GUST_OG.txt","TUSON_OG_reported_by_2020plus.txt")

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

methods<-c("2020_OG.txt","TUSON_OG.txt")

for(i in methods){

		prediction<-read.table(paste("topgenes/",i,sep=""),header=F);
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

result<-cbind(c("DORGE","20/20+","GUST","TUSON","20/20+*","TUSON*"),result)

result<-as.data.frame(result)
colnames(result)<-c("Method","#","Sn","Sp","Precision","Accuracy")
result$Accuracy<-as.numeric(result$Accuracy)
result<-result[c(order(-result[1:4,]$Accuracy),c(5,6)),]

write.table(result,"OG_gene_performance.txt",sep="\t",quote=F,row.names=F,col.names=T)