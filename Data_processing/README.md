# Data processing used in the DORGE paper
![image](https://github.com/biocq/DORGE/blob/master/DORGE_logo.svg)


## Availability
Codes and data used for processing data as well as outputs are available at Figshare website. See specific folders for details.

## Dependency

R (https://cran.r-project.org/)

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.


## Illustration
The illustration of the data processing is presented in README.md file in individual folders, including:

*  Gene_annotation/README.md
*  Epigenetics_processing/README.md
*  Genomics_processing/README.md
*  Mutation_processing/README.md
*  Phenotype_processing/README.md
*  Features/README.md


The features can by combined by using github_feature_combination.java provided here.

## Instructions

**1. Download datasets and place in the proper hierarchy.**

The layout of data folders that we expect:

|- **Data_processing** (root folder)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Gene_annotation (Gene annotation based on Ensembl, GENCODE, COSMIC and other resources)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Epigenetics_processing (Pipeline for processing epigenetic datasets)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Genomics_processing (Pipeline for processing genomic datasets)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Mutation_processing (Pipeline for processing mutation datasets)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Phenotype_processing (Pipeline for processing phenotype datasets)
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Features (Including a script for generating training feature profile)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- github_feature_combination.java

**2. `cd` to individual folders under DORGE_pipeline and follow the instructions (README.md) to run the codes for specific feature categories.**

**3. Script for combining features**
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

Output: All_features_training.csv (Training feature profile including 75 features used for model training (https://github.com/biocq/DORGE))

## Changelog
*  March 25. Update the scripts and fix bugs.
*  April 27. Update the README page.
