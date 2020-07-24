# Codes used for Gene annotation


## Codes to annotate genes can be found at:

https://figshare.com/s/00e9fd6341023a6b8770

## Dependency

Java v1.8.0

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

## Steps

1. Ensembl synonyms to gene symbol
```
javac step1_ensembl_synonyms_to_symbol.java
java step1_ensembl_synonyms_to_symbol
```

  Input: external_synonym.txt #downloaded from ensembl-annotation
  
  Input: xref.txt #Downloaded from The Ensembl Xref System
  
  Output: ensembl_synonymous_to_symbol.txt


2. Genes annotated by COSMIC annotation
```
javac step2_gene_merge.java
java step2_gene_merge
```
Input: H3K4me3_width.all_in_one.xls
Input: TUSON_all_genelist.txt
Output: compiled_genelist.txt #(Gene list compiled from TUSON and Kaifu Chen's Paper)

3. Gene name list generation
```
javac step3_gene_name_generation.java
java step3_gene_name_generation
```
  Input: Gene_set_new.txt

  Input: TUSON_input_new_sorted.txt # From COSMIC gene list
  
  Output: compiled_genelist_genename.txt # (Gene symbol list consistent with Gene_set_new.txt)
  
  Output: gene_set_sorted.txt #(Gene symbol list consistent with mutation with same sequence)

4. Gene promoter and genebody BED file generation
```
javac Step4_gene_promoter_genebody_BED_generation.java
java Step4_gene_promoter_genebody_BED_generation
```
  Input: compiled_genelist_genename.txt #(Gene symbol list consistent with Gene_set_new.txt)
  
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
  
  Input: gencode.v25lift37.annotation.gtf #(GENCODE annotation)
  
  Input: gencode.v19.annotation.gtf #(GENCODE annotation)

  Output: compiled_genelist_promoter_hg19.bed #(Gene promoters defined as [-500, 1000] around TSSs)

  Output: compiled_genelist_promoter_250_hg19.bed #(Gene promoters defined as [-500, 250] around TSSs)

  Output: compiled_genelist_genebody_hg19.bed #(Gene promoters defined as 500 downstream of TSSs towards to TTSs)


5. Gene full-length BED file generation
```
javac Step5_gene_fulllength_BED_generation.java
java Step5_gene_fulllength_BED_generation
```

  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Input: compiled_genelist_genename.txt #(Gene symbol list consistent with Gene_set_new.txt)
  
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
  
  Input: gencode.v25lift37.annotation.gtf #(Gencode v25 hg19 GTF annotation)
  
  Input: gencode.v19.annotation.gtf #(Gencode v19 hg19 GTF annotation)

  Output: compiled_gene_pos_hg19.bed

6. Gene name and alias mapping file generation
```
javac Step6_COSMIC_to_Gene_Symbol.java
java Step6_COSMIC_to_Gene_Symbol
```
  Input: compiled_genelist_genename.txt #(Gene symbol list consistent with Gene_set_new.txt)
  
  Input: CosmicTranscripts.tsv #(Transcript anntation file downloaded from COSMIC website)
  
  Input: Ensembl_TransID_v87_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)
  
  Input: Ensembl_TransID_v92_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)
  
  Output: Cosmic2Gene_name.txt #(Gene name and alias mapping file)