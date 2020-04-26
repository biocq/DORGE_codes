## Codes used in the DORGE paper
![image](https://github.com/biocq/DORGE/blob/master/www/DORGE_logo.svg)
### Codes for reproducing the results of the paper
Codes and necessary data used in the DORGE paper for generating figures and tables are available at this website.

### Dependency

R (https://cran.r-project.org/)

Java (Optional if only you are only interested in regenerating figures)

Dependency on R packages can be found in R files and can be installed automatically.

Sun Java SE Development Kit 7 or higher (https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is also needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

### Usages

Take Fig1_AB_code.R in folder Figure_1 for an example.

Option 1. In terminal/console

Change working directory to the folder containing Fig1_AB_code.R.

run "Rscript Fig1_AB_code.R".

Option 2. Run in R GUI

Run setwd("directory to Figure_1") to change to Figure_1 folder. You can use source("Fig1_AB_code.R", echo = TRUE) to run the entire script.

Option 3. Run in Rstudio

Under session> Set Working directory> Choose directory... to change to Figure_1 folder or setwd("directory to Figure_1"). You can use source("Fig1_AB_code.R", echo = TRUE) to run the entire script.

Option 4. Batch mode

We also provide a convenient executable script (script_batch) for running R codes in a batch mode. Please set permission (i.e. chmod a+x script_batch) before excecuting it.

### Layout of the website

	 |-Figure_1
   |---Raw_figures                        Unprocessed raw figures related to Figure 1 and Figure S1
   |---data					                      Processed data used for figure generation
   |---Figure_1.pdf                       Assembled Figure 1
   |---Figure S1 related to Figure 1.pdf  Assembled Figure S1
   |-Figure_2
   |---Raw_figures                        Unprocessed raw figures related to Figure 2 and Figure S2
   |---data					                      Processed data used for figure generation
   |---Figure_2.pdf                       Assembled Figure 2
   |---Figure S2 related to Figure 2.pdf  Assembled Figure S2
   |-Figure_3
   |---Raw_figures                        Unprocessed raw figures related to Figure 3 and Figure S3
   |---data					                      Processed data used for figure generation
   |---Figure_3.pdf                       Assembled Figure 3
   |---Figure S3 related to Figure 3.pdf  Assembled Figure S3
   |-Figure_4
   |---Raw_figures                        Unprocessed raw figures related to Figure 4 and Figure S4
   |---data					                      Processed data used for figure generation
   |---Figure_4.pdf                       Assembled Figure 4
   |---Figure S4 related to Figure 4.pdf  Assembled Figure S4
   |-Figure_5
   |---Raw_figures                        Unprocessed raw figures related to Figure 5 and Figure S5-S7
   |---data					                      Processed data used for figure generation
   |---Figure_5.pdf                       Assembled Figure 5
   |---Figure S5 related to Figure 5.pdf  Assembled Figure S5
   |---Figure S6 related to Figure 5.pdf  Rearranged Figure S6
   |---Figure S7 related to Figure 5.pdf  Rearranged Figure S7
   |-Tables
   |---Data_file_S1                       Unprocessed raw tables that are used to make the final Data file S1
   |---Data_file_S2                       Unprocessed raw tables that are used to make the final Data file S2
   |---Table_1                            Unprocessed raw table of Table 1
   script_batch				                    A convenient executable script (script_batch) for running R codes in a batch mode
   Gene_set_new.txt		                    Gene annotation file through the paper
   DORGE_prediction.txt                   DORGE prediction results
   All_features.csv                       The feature profile

## Citation

Please cite the paper (DORGE: Discovery of Oncogenes and Tumor SuppressoR Genes Using Genetic and Epigenetic Features, now submitted to Science Advances) if the resources are used elsewhere.


## Changelog

*  April 4. Update the README page.

*  April 25. Update the codes.

## Issues

Codes have been tested, if you encounter errors, please check the dependency or restart the R environment, given methods in packages may conflict and unexpected errors may occur. If you still encounter any problems, please [file an issue](https://github.com/biocq/DORGE_paper/issues) along with a detailed description.