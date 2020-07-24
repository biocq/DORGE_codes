# Codes used for Phenotype feature calculation


  
## Codes to process Phenotype data can be found at:

https://figshare.com/s/da23c6485dfde9632f73

## Dependency

Java v1.8.0

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

## Steps


1. Gene expression data analysis.

```
cd Gene_expression
gunzip CosmicCompleteGeneExpression.tsv.gz
gunzip CosmicSample.tsv.gz

javac COSMIC_gene_expression_data_tissues.java
java COSMIC_gene_expression_data_tissues

javac -cp ../../Genomics_processing/replication_time/commons-math3-3.6.1.jar COSMIC_gene_expression_sample_statistics.java
java -Xmx2048m -cp ../../Genomics_processing/replication_time/commons-math3-3.6.1.jar:. COSMIC_gene_expression_sample_statistics

javac Gene_expression_value_calculation.java
java Gene_expression_value_calculation
cd ..
```
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)

  Input: gencode.v19.annotation.gtf #(Gencode v19 GTF annotation)

  Input: CosmicCompleteGeneExpression.tsv #(COSMIC Gene expression data v90)
  
  Internal: Tissues_expression_data.txt #(Samples in gene expression data)

  Output: gene_expression_quant.txt #(Mean expression levels)

2. CRISPR-screening data analysis (Time-consuming, ~ 15 min).
```
cd CRISPR_screen
javac -cp ../../Genomics_processing/replication_time/commons-math3-3.6.1.jar avana_median_dependency_score_calculation.java
java -Xmx1524m  -cp ../../Genomics_processing/replication_time/commons-math3-3.6.1.jar:. avana_median_dependency_score_calculation
cd ..
```
  Input: Achilles_gene_effect.csv #(CRISPR-screening dependency data (DepMap Public 20Q1), downloaded from https://depmap.org/portal/download/)

  Internal: Achilles_gene_effect_transposed.csv #(Transposed dependency matrix and fill NA with the median dependency value in specific cell lines)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)

  Output: CRISPR_screening_dependency_scores.txt #(Mean Independency and dependency scores)


3. VEST score.
```
cd VEST
javac VEST_feature_score_calculation.java
java VEST_feature_score_calculation
cd ..
```
  Input:  VEST_raw.txt #(from step 9b result in mutation processing)

  Output: VEST_values.txt #(Mean VEST scores)
