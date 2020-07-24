# Pipeline used to process data in the DORGE paper
![image](https://github.com/biocq/DORGE/blob/master/DORGE_logo.svg)


## Availability
Codes and data used for processing data as well as outputs are available at Figshare website. See specific folders for details.

## Dependency

R (https://cran.r-project.org/)

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.


## Illustration
The illustration of the pipeline is presented at individual folders, including:

*  Gene_annotation (Gene annotation based on Ensembl, GENCODE, COSMIC and other resources)
*  Epigenetics_processing (Pipeline for processing epigenetic datasets)
*  Function_processing (Pipeline for processing function datasets)
*  Genomics_processing (Pipeline for processing genomic datasets)
*  Mutation_processing (Pipeline for processing mutation datasets)
*  Phenotype_processing (Pipeline for processing phenotype datasets)
*  Features (Including a script for generating training feature profile)


The features can by combined by using github_feature_combination.java provided in this website.

## Instructions

**1. Download datasets and place in the proper hierarchy.**

The layout of data folders is expected:

|- **DORGE_pipeline** (folder containing codes and data for generating feature profiles used by the DORGE paper for model training and figure/tables generation)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Gene_annotation (Gene annotation based on Ensembl, GENCODE, COSMIC and other resources)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Epigenetics_processing (Pipeline for processing epigenetic datasets)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Function_processing (Pipeline for processing function datasets)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Genomics_processing (Pipeline for processing genomic datasets)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Mutation_processing (Pipeline for processing mutation datasets)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Phenotype_processing (Pipeline for processing phenotype datasets)
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Features (Including a script for generating training feature profile)

&nbsp;&nbsp;&nbsp;&nbsp;|--- Features (Including a script for generating training feature profile)

&nbsp;&nbsp;&nbsp;&nbsp;|- github_feature_combination.java

**2. Go to individual folders under DORGE_pipeline and follow the instructions (README.md) to run the codes for specific feature categories.**

**3. Code for combining features**
```
javac github_feature_combination.java
java github_feature_combination
```

Output: All_features.csv (Complete feature profile including 97 features in the Features folder)

**4. Code for producing training feature subset**
```
cd Features
Rscript feature_subset_generation.R
cd ..
```

Output: All_features_training.csv (Training feature profile including 75 features in the Features folder)

## Citation

Please cite the paper (DORGE: Discovery of Oncogenes and Tumor SuppressoR Genes Using Genetic and Epigenetic Features, https://www.biorxiv.org/content/10.1101/2020.07.21.213702v1) if the resources are used elsewhere.


## Changelog
*  March 25. Update the scripts and fix bugs.
*  April 27. Update the README page.

## Issues

If you encounter any problems, please [file an issue](https://github.com/biocq/DORGE_pipeline/issues) along with a detailed description.
