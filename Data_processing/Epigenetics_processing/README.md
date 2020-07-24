# Codes used for Epigenetics feature calculation

## Codes to process evaluation data can be found at:

https://figshare.com/s/8b8fa54f490344c31620

## Dependency

Java v1.8.0, bedtools (https://bedtools.readthedocs.io/en/latest/), Subread (http://subread.sourceforge.net/)

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

## Steps

1. COSMIC methylation data analysis (including promoter and genebody methylation)
```
cd methylation
gunzip CosmicCompleteDifferentialMethylation.tsv.gz
javac genebody_promoter_methylation_processing.java
java genebody_promoter_methylation_processing
cd ..
```
  Input: genebody_450k_probe_overlap.txt #(Gene body to probe mapping)
  
  Input: promoter_450k_probe_overlap.txt #(Gene promoter to probe mapping)
  
  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
  
  Input: CosmicCompleteDifferentialMethylation.tsv # (COSMIC 450K methylation data v90)
  
  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)

  Output: genebody_differential_methyl_COSMIC.txt # (Genebody methylation information in cancer and normal samples)

  Output: promoter_differential_methyl_COSMIC.txt # (Promoter methylation information in cancer and normal samples)

2. Gene-body canyon analysis

2. (a) Generation of canyon gene BED file
```
cd genebody_canyon
bedtools intersect -wo -a Canyon_Hyper.txt -b ../../Gene_annotation/Step4_gene_promoter_genebody_BED_generation/compiled_genelist_genebody_hg19.bed > Canyon_with_genebody.bed
javac genebody_canyon_processing.java
java genebody_canyon_processing
bedtools sort -i Canyon_with_genebody_processed.bed > Canyon_with_genebody_processed.sorted.bed
```
 
2. (b) Association of genebody canyon genes with 450K methylation information
```
javac WGBS_canyon_processing.java
java WGBS_canyon_processing
cd ..
```
  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Internal: Canyon_with_genebody_processed.bed (Internal result)

  Output: Canyon_methylation.txt (450K methylation information for genebody canyon genes)

3. Histone modification data analysis

3. (a) Histone modification data preparation
```
cd codes
unzip H3K27ac.zip #uncompress internal results
unzip H3K27me3.zip
unzip H3K36me3.zip
unzip H3K4me1.zip
unzip H3K4me2.zip
unzip H3K4me3.zip
unzip H3K79me2.zip
unzip H3K9ac.zip
unzip H3K9me2.zip
unzip H3K9me3.zip
unzip H4K20me1.zip
chmod a+x all_code.txt # Very time-consuming; Download and process histone modification data
./all_code.txt
cd ..
```
   Input: all_code.txt
   Output: histone modification processed files within codes folder
3. (b) Histone modification quantification
```
cd histone_modification
javac histone_modification_data_processing.java
java histone_modification_data_processing
cd ..
```
  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Input: compiled_gene_pos_hg19.bed #(From step 5 from Gene_annotation folder)
  
  Input: histone modification processed files within codes folder

  Output: ENCODE_H3K4me1_peak_info.txt #(Mean peak length, Percentage of broad peaks and Height of peaks)

  Output: ENCODE_H3K4me2_peak_info.txt #(Mean peak length, Percentage of broad peaks and Height of peaks)

  Output: ENCODE_H3K4me3_peak_info.txt #(Mean peak length, Percentage of broad peaks and Height of peaks)

  Output: ENCODE_H3K9ac_peak_info.txt #(Mean peak length, Percentage of broad peaks and	Height of peaks)

  Output: ENCODE_H3K9me2_peak_info.txt #(Mean peak length, Percentage of broad peaks and Height of peaks)

  Output: ENCODE_H3K9me3_peak_info.txt #(Mean peak length, Percentage of broad peaks and Height of peaks)

  Output: ENCODE_H3K27ac_peak_info.txt #(Mean peak length, Percentage of broad peaks and Height of peaks)

  Output: ENCODE_H3K27me3_peak_info.txt #(Mean peak length, Percentage of broad peaks and	Height of peaks)

  Output: ENCODE_H3K36me3_peak_info.txt #(Mean peak length, Percentage of broad peaks and	Height of peaks)

  Output: ENCODE_H3K79me2_peak_info.txt #(Mean peak length, Percentage of broad peaks and	Height of peaks)

  Output: ENCODE_H4K20me1_peak_info.txt #(Mean peak length, Percentage of broad peaks and	Height of peaks)
  
  
4. Super enhancer.
```
cd Super_enhancer
unzip all_hg19_bed.zip
chmod a+x code_1_sort.txt
./code_1_sort.txt #(sort BED files)
chmod a+x code_2_bedtools.txt
./code_2_bedtools.txt #(Link super enhancers to genes)
javac super_enhancer_result_processing.java
java super_enhancer_result_processing
cd ..
```
  Input:  BED files in the all_hg19_bed #(Downloaded from dbSUPER http://bioinfo.au.tsinghua.edu.cn/dbsuper/)

  Input:  compiled_gene_pos_hg19.sorted.bed #(Gene position annotation hg19 version)

  Internal: Super_enhancer_overlap.txt #(For diagnosis)

  Output: super_enhancer_stats.txt #(Cell-lines that are associated with super enhancers for specific genes)


5. Cell cyle replication timing data

```
cd replication_time
./script # Replication timing quantification (Very time-consuming)
```
  Input: Repli-seq Raw files
  
  Input: compiled_gene_pos_hg19.bed #(Gene position annotation hg19 version)
  
  Output: quant folder (Replication timing quantification files)

```
javac -cp commons-math3-3.6.1.jar replication_timing.java 
java -cp commons-math3-3.6.1.jar:. replication_timing
cd ..
```
  Input: Replication timing quantification files (quant folder)
  
  Input: Gene_set_new.txt #(Gene annotation information used by this study)
  
  Output: Replication_timing_percent.txt (Replication timing S50 score profile)
  
  Output: Replication_timing_stats.txt
  