# Codes used in the DORGE paper
![image](https://github.com/biocq/DORGE/blob/master/DORGE_logo.svg)
## Codes for reproducing the results of the paper
Codes and necessary data used in the DORGE paper for generating figures and tables are available at this website.

## Dependency

*  **R** 4.0.0 or higher (https://cran.r-project.org/)

Dependency on R packages can be found in R files and can be installed automatically. Versions of R packages available for troubleshooting are shown below.

R Package	Versions  
**circlize**	0.4.9  
**ComplexHeatmap**	2.4.1  
**cowplot**	1.0.0  
**data.table**	1.12.8  
**dplyr**	0.8.5  
**ggplot2**	3.3.0  
**ggpubr**	0.2.5  
**ggrepel**	0.8.2  
**ggsignif**	0.6.0  
**Hmisc**	4.4-0  
**plyr**	1.8.6  
**PRROC**	1.3.1  
**reshape2**	1.4.4  
**scales**	1.1.0  
**VennDiagram**	1.6.20  

*  **Cytoscape** 3.8.0 or higher (Optional; If you are only interested in viewing network session file of Cytoscape)

*  **Java SDK** 8 or higher (https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html)

Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

## Usages

### Figures/tables generation

Take Fig1_AB_code.R in folder Figure_1 for an example.

Option 1. In terminal/console

Change working directory to the folder containing Fig1_AB_code.R.

run "Rscript Fig1_AB_code.R".

Option 2. Run in R GUI

Run setwd("directory to Figure_1") to change to Figure_1 folder. You can use source("Fig1_AB_code.R", echo = TRUE) to run the entire script.

Option 3. Run in Rstudio

Under session> Set Working directory> Choose directory... to change to Figure_1 folder or setwd("directory to Figure_1"). You can use source("Fig1_AB_code.R", echo = TRUE) to run the entire script.

Option 4. Batch mode

We also provide a convenient executable script (**script_batch**) for running R codes in a batch mode. Please set permission (i.e. chmod a+x script_batch) before executing it.

### Processing data to generate feature profile

Instructions are included in Data_processing folder

## Layout of the website
<br/>
|-Figure_1<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Raw_figures                        			Unprocessed raw figures related to Figure 1 and Figure S1<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---data			                      		Processed data used for figure generation<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure_1.pdf                       			Assembled Figure 1<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure S1 related to Figure 1.pdf  			Assembled Figure S1<br/>
|-Figure_2<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Raw_figures                        			Unprocessed raw figures related to Figure 2 and Figure S3<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---data			                      		Processed data used for figure generation<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure_2.pdf                       			Assembled Figure 2<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure S3 related to Figure 2.pdf  			Assembled Figure S3<br/>
|-Figure_3<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Raw_figures                        			Unprocessed raw figures related to Figure 3 and Figure S4<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---data			                      		Processed data used for figure generation<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure_3.pdf                       			Assembled Figure 3<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure S4 related to Figure 3.pdf  			Assembled Figure S4<br/>
|-Figure_4<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Raw_figures                        			Unprocessed raw figures related to Figure 4 and Figure S5-S6<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---data			                      		Processed data used for figure generation<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure_4.pdf                       			Assembled Figure 4<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure S5 related to Figure 4.pdf  			Assembled Figure S5<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Figure S6 related to Figure 4.pdf  			Assembled Figure S6<br/>
|-Tables<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Data_file_S1                       			Raw tables that are used to make the final Data file S1<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Data_file_S2                       			Raw tables that are used to make the final Data file S2<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Table_1                            			Raw table of Table 1<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Data file S1.xlsx                  			The final Data file S1<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Data file S2.xlsx                  			The final Data file S2<br/>
|-Data_processing<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Epigenetics_processing             			Instruction of epigenetic data processing<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Features		            			 			Included combined feature profiles and codes to filter features<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Gene_annotation             			 			Instruction of gene annotation processing<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Genomics_processing             	 			Instruction of genomics data processing<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Mutation_processing             	 			Instruction of mutational data processing<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---Phenotype_processing             	 			Instruction of phenotype data processing<br/>
&nbsp;&nbsp;&nbsp;&nbsp;|---github_feature_combination.java    			Java program to combine features<br/>

script_batch		                    								A convenient executable script (script_batch) for running R codes in a batch mode<br/>
Gene_set_new.txt		                    								Gene annotation file through the paper<br/>
DORGE_prediction.txt                    								DORGE prediction results<br/>
All_features.csv                        								The feature profile<br/>

## Citation

Please cite the paper (DORGE: Discovery of Oncogenes and Tumor SuppressoR Genes Using Genetic and Epigenetic Features, https://www.biorxiv.org/content/10.1101/2020.07.21.213702v1) if the resources are used elsewhere.


## Changelog
*  April 4. Update the README page.
*  April 25. Update the codes.
*  May 7. Update the codes.
*  July 20. Update the codes.

## Issues

Codes have been tested in R 4.0.0 and Java 1.8.0. If you encounter errors, please check the dependency or restart the R environment, given packages may change in future and methods in packages may conflict and unexpected errors may occur. If you still encounter any problems, please [file an issue](https://github.com/biocq/DORGE_paper/issues) along with a detailed description.
