###########################################
###### DGB data processing  ###############
###########################################

###### Generate tables of expression data for DORGE-predicted TSG/OGs that are up-/down-regulated by drug perturbations from the Connectivity Map (CMap) collection in the Drug Gene Budger (DGB) web portal.

###### Input: ../../Figure_5/data/DGB/DGB_creeds_downregulated_combined.txt: Processed expression data for DORGE-predicted OGs that are down-regulated by drug perturbations from the CRowd Extracted Expression of Differential Signatures (CREEDS) collection in the Drug Gene Budger (DGB) web portal.
###### Input: ../../Figure_5/data/DGB/DGB_creeds_upregulated_combined.txt: Processed expression data for DORGE-predicted TSGs that are up-regulated by drug perturbations from the CRowd Extracted Expression of Differential Signatures (CREEDS) collection in the Drug Gene Budger (DGB) web portal.

###### Output: : DGB_creeds_downregulated_formatted.txt: Formatted expression data for DORGE-predicted OGs that are down-regulated by drug perturbations from the CRowd Extracted Expression of Differential Signatures (CREEDS) collection in the Drug Gene Budger (DGB) web portal.
###### Output: : DGB_creeds_upregulated_formatted.txt: Formatted expression data for DORGE-predicted TSGs that are up-regulated by drug perturbations from the CRowd Extracted Expression of Differential Signatures (CREEDS) collection in the Drug Gene Budger (DGB) web portal.

options(warn=-1)
################################### Data file S1: Formatted gene expression UP/down data from CREEDS collection in the Drug Gene Budger (DGB) web portal
 ###################################
#setwd("/Users/jlyu/Box\ Sync/TSGOG_Project/SA_sub/github/DORGE_paper/DORGE_codes/Tables/Data_file_S1/");

data <- read.table("../../Figure_5/data/DGB/DGB_creeds_downregulated_combined.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(data)<-c("Gene","Gene ID","Gene type","DORGE_OG_score","Drug Name","CREEDS ID","GEO ID","DrugBank ID","P-value","Q-value","Log2(fold Change)","Specificity");
write.table(data[,c("Gene","Gene ID","Gene type","Drug Name","CREEDS ID","GEO ID","DrugBank ID","P-value","Q-value","Log2(fold Change)")], "DGB_creeds_downregulated_formatted.txt",sep="\t",quote=F,row.names=F,col.names=T)

data <- read.table("../../Figure_5/data/DGB/DGB_creeds_upregulated_combined.txt", header=T, sep="\t",fill=TRUE,quote = "")
colnames(data)<-c("Gene","Gene ID","Gene type","DORGE_TSG_score","Drug Name","CREEDS ID","GEO ID","DrugBank ID","P-value","Q-value","Log2(fold Change)","Specificity");
write.table(data[,c("Gene","Gene ID","Gene type","Drug Name","CREEDS ID","GEO ID","DrugBank ID","P-value","Q-value","Log2(fold Change)")], "DGB_creeds_upregulated_formatted.txt",sep="\t",quote=F,row.names=F,col.names=T)
