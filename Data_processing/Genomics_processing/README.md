# Codes used for Genomics feature calculation

## Codes to process Genomics data can be found at:

https://figshare.com/s/2e5c4a59673385545f6f

## Dependency

Java v1.8.0

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

## Steps
  
1. NCBOOST_features feature processing

```
cd NCBOOST_features
javac NCBoost_geneDB.java
java NCBoost_geneDB
cd ..
```
  Input: NCBoost_geneDB.tsv
  
  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
  
  Output: NCBoost_geneDB_processed.txt #(Several genomic and evoluation related features)
  
2a. CNA feature processing

```
cd CNA
gunzip CosmicCompleteCNA.tsv.gz #v90
javac CNA_COSMIC_processing.java
java  -Xmx1024m CNA_COSMIC_processing #(Time-consuming; ~ 30min)
```
  Input: CosmicCompleteCNA.tsv #(COSMIC CNA raw data v90)
  
  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
  
  Output: CNA_gene_CNA_status.txt #(Gain or loss information)
  
  Output: CNA_gene_sampleNum.txt #(Total sample numbers associated with specific genes)

2b. CNA feature calculation
```
javac CNA_feature_calculation.java
java CNA_feature_calculation
cd ..
```
  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
  
  Input: CNA_gene_CNA_status.txt #(From last step)
  
  Input: CNA_gene_sampleNum.txt #(From last step)
  
  Output: CNA_freq.txt #(Amplification and Deletion percentage)
