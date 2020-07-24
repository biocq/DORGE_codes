## Codes used for mutation feature calculation

## Codes to process evaluation data can be found at:

https://figshare.com/s/fb7b1f4db437f3eb1536

## Dependency

Java v1.8.0

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

## Steps

1. COSMIC mutation merging
```
cd step1_COSMIC_Somatic_mutation_filter
unzip CosmicMutantExport.tsv.zip #(Raw data, compressed to upload)
unzip Mutation_COSMIC_2020.txt.zip #(Result, compressed to upload)
javac step1_COSMIC_Somatic_mutation_filter.java
java step1_COSMIC_Somatic_mutation_filter
cd ..
```
  Input: Cosmic Mutation files #(CosmicMutantExport.tsv.zip downloaded from COSMIC)
 
  Output: Mutation_COSMIC_2020.txt #(Merged COSMIC mutation MAF file)

2. TCGA mutation merging
```
cd step2_TCGA_Somatic_mutation_merge
cd muse
unzip SNV_1.zip #(Raw data, compressed to upload)
unzip SNV_2.zip #(Raw data, compressed to upload)
cd ..
unzip Mutation_TCGA_2020.txt.zip #(Result, compressed to upload)
javac step2_TCGA_Somatic_mutation_merge.java
java step2_TCGA_Somatic_mutation_merge
cd ..
```
  Input: TCGA Mutation files #(TCGA pan-cancer TCGA.*.muse maf files in the muse folder)
 
  Output: Mutation_TCGA_2020.txt #(Merged COSMIC mutation file in 20/20+ format)


3. COSMIC and TCGA mutation merging
```
cd step3_TCGA_COSMIC_Somatic_mutation_merge
unzip Mutation_merged_2020.txt.zip #(Result, compressed to upload)
unzip Mutation_merged.txt.zip #(Result, compressed to upload)
javac step3_TCGA_COSMIC_Somatic_mutation_merge.java
java step3_TCGA_COSMIC_Somatic_mutation_merge
cd ..
```
  Input: Mutation_COSMIC_2020.txt #(from step 1)

  Input: Mutation_TCGA_2020.txt  #(from step 2)
 
  Output: Mutation_merged_2020.txt #(Merged and nonredundant mutation file in 20/20+ format)
 
  Output: Mutation_merged.txt #(Merged and nonredundant mutation file in TUSON format)

4. BED file generation #(both hg19 and Grch38, for 20/20+ software)
```
cd step4_bed12_file_generation
unzip Mutation_merged_2020_refined_hg19_final.txt.zip #(Result, compressed to upload)
unzip Mutation_merged_2020_refined_hg19.txt.zip #(Result, compressed to upload)
javac step4_bed12_file_generation.java
java step4_bed12_file_generation
cd ..
```
  Input: various Ensembl gene/transcript ID annotation files #(Downloaded or processed from Ensembl raw files)

  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)

  Input: CosmicTranscripts.tsv #(Transcript anntation file downloaded from COSMIC website)

  Output: genename_hg19_for_2020.bed12 #(Generate gene bed file for genes mapped to hg19 genome version)
 
  Output: genename_Grch38_for_2020.bed12 #(Generate gene bed file for genes that can only been mapped to Grch38 genome version)

  Input: genename_hg19_for_2020_comverted_from_Grch38_with_fails.bed #(CrossMap.py unmapped genes #(from Step 4), used in generating formated Mutation file)

  Input: Mutation_merged_2020_refined_hg19.txt
 
  Output: Mutation_merged_2020_refined_hg19_final.txt #(formatted Mutation file)

5. Convert genename_Grch38_for_2020.bed12 from hg38 to hg19
```
cd step5_genome_version_conversion
CrossMap.py bed hg38ToHg19.over.chain.gz ../step4_bed12_file_generation/genename_Grch38_for_2020.bed12 > genename_hg19_for_2020_comverted_from_Grch38.bed
cd ..
```
  CrossMap can be downloaded from http://crossmap.sourceforge.net/

  Input: hg38ToHg19.over.chain.gz (Chain files from hg38 (GRCh38) to hg19)
 
  Output: genename_hg19_for_2020_comverted_from_Grch38.bed

6. BED file combination
```
cd step6_combine_bed12
javac step6_combine_bed12.java
java step6_combine_bed12
cd ..
```
  Input: genename_hg19_for_2020.bed12 #(from step 4)

  Input: genename_hg19_for_2020_comverted_from_Grch38.bed #(from step 5)

  Input: Mutation_merged_2020_refined_hg19.txt #(from step 4)
 
  Output: final_genes_hg19_2020.bed12

7. Mutation files prepared as input of 20/20+, CRAVAT and PolyPhen-2
```
cd step7_mutation_CRAVAT
unzip gencode.v20.annotation.gtf.zip
unzip Mutation_merged_2020_refined.txt.zip #(Result, compressed to upload)
unzip Mutation_merged_2020_refined_hg19.txt.zip #(Result, compressed to upload)
javac step7_mutation_CRAVAT.java
java step7_mutation_CRAVAT
cd ..
```
  Input: Mutation_merged_2020.txt (Output from step 3)

  Input: gencode_genename.v28.hg19.bed

  Input: mutation_base_from_coordinate.txt #(mutation base from coordinate generated from Step 8)

  Input: gencode.v20.annotation.gtf #(Gencode V20 annotation)

  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)

  Input: strand_need_correction.txt #(Some mutation record has wrong chain information, refused by online website, need correcting)

  Input: Mutation_hg19.bed #(Mutation genome position hg38 to hg19 conversion information from LiftOver)
 
  Output: Mutation_merged_CRAVAT.txt #(Mutation file GRCh38 for CRAVAT webserver)
 
  Output: Mutation_merged_2020_refined.txt #(Mutation file for 20/20+ software)
 
  Output: Mutation_merged_2020_refined_hg19.txt #(Mutation file for 20/20+)
 
  Output: Mutation_Gcrh38.bed #(Mutation file with Grch38 genome version)
 
  Output: Mutation_hg19_input_polyphen.txt #(Mutation file for PolyPhen-2 webserver)

8. Sort three mutation files and convert a BED12 file to hg19

```
cd step8_mutation_sort
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
unzip hg38.fa.gz
bedtools getfasta -fi hg38.fa -bed Mutation_Gcrh38.bed -fo ../step7_mutation_CRAVAT/mutation_base_from_coordinate.txt -tab -s
sort -k2,2 ../step7_mutation_CRAVAT/Mutation_merged_CRAVAT.txt > Mutation_merged_CRAVAT.sorted.txt
sort -k1,1n ../step7_mutation_CRAVAT/Mutation_merged_2020_refined.txt > Mutation_merged_2020_refined.sorted.txt
sort -k1,1n ../step7_mutation_CRAVAT/Mutation_merged_2020_refined_hg19.txt > Mutation_merged_2020_refined_hg19.sorted.txt
cd ..
```
9. (a) CRAVAT VEST score generation
Mutation_merged_CRAVAT.txt submitted to CRAVAT website (http://www.cravat.us/CRAVAT/)

9. (b) CRAVAT VEST score merge
```
cd step9_VEST_feature_combination
unzip VEST_results.zip
javac step9_VEST_feature_combination.java
java step9_VEST_feature_combination
cd ..
```
  Input: CRAVAT_results folder generated by the CRAVAT website (http://www.cravat.us/CRAVAT/)

  Input: Gene_set_new.txt #(Gene set annotation)

  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
 
  Output: VEST_raw.txt #(VEST scores for all mutations)

10. (a) Mutation_hg19_input_polyphen.txt (from step 7) submitted to PolyPhen-2 web server

10. (b) PolyPhen-2 input file preparation, rerun by fixing wrong mutation information based on PolyPhen-2 logs and integrate PolyPhen-2 scores based on two round PolyPhen-2 running (polyphen_input2.txt and polyphen_input3.txt).
```
cd step10_polyphen2_input_preparation
unzip *.zip
javac step10_polyphen2_input_preparation.java
java step10_polyphen2_input_preparation
cd ..
```
  Input: Mutation_hg19_input_polyphen.txt; 

  Input: Mutation_hg19_input_polyphen_old.txt

  Input: pph2-log_new.txt pph2-log_1.txt #(PolyPhen-2 error information)
 
  Output: polyphen_input3.txt #(3rd round submission)
 
  Output: polyphen_input2.txt #(2nd round submission)

  Output: pph2-short_all.txt #(PolyPhen-2 results)
 
  Output: pph2-short_new.txt #(PolyPhen-2 results)

`cat pph2-short_all.txt pph2-short_new.txt > pph2-short_all_new.txt` #(Final PolyPhen-2 results by Combining 1st round PolyPhen-2 results with 2nd round results based on corrected mutation annotation)

11. (a) CRAVAT mutation files genome conversion, input to CRAVAT server at http://hg19.cravat.us/CRAVAT/
```
cd step11_CRAVAT_to_hg19
javac CRAVAT_to_hg19.java
java CRAVAT_to_hg19
cd ..
```
  Input: Mutation_hg19.bed #(from Step 7 output; genome mapping file)

  Input: Mutation_merged_CRAVAT.sorted.txt #(from Step 8 output)
 
  Output: Mutation_merged_CRAVAT_hg19.txt

11. (b) Submit Mutation_merged_CRAVAT_hg19.txt to CRAVAT
CRAVAT website (http://hg19.cravat.us/CRAVAT/ with SnvGet options)


12. TUSON input preparation
```
unzip pph2-short_all_processed.txt.zip
unzip newest_Mutation_merged_withoutPolyphen_scores.txt.zip
cd step12_TUSON_input_prepare
javac step12_TUSON_input_prepare.java
java step12_TUSON_input_prepare
cd ..
```
  Input: pph2-short_all_new.txt #(from Step 10 output)

  Input: Mutation_hg19.bed #(from Step 7 output; genome mapping file)

  Input: Gene_set_new.txt

  Input: gene_set_sorted.txt

  Input: TUSON_input_new_sorted.txt (from output of Step3_gene_name_generation)
 
  Output: gene_set_sorted.txt

  Input: gencode.V20_ENSEMBL_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)

  Input: Ensembl_v87_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)

  Input: Ensembl_v92_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)

  Input: Ensembl_TransID_v87_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)

  Input: Ensembl_TransID_v92_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)

  Input: CosmicTranscripts.tsv #(Transcript anntation file downloaded from COSMIC website)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)

  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)

  Input: Mutation_merged_2020_refined.txt #(from Step 6 output)

  Input: all_NA.txt #(imputation for some genes that are all NAs in pph2 score)
 
  Output: TUSON_input_new_revised.txt #(Final TUSON input file)
  
  Internal: all_pph.txt #(Internal result)

  Output: TUSON_input_new.txt #(Internal result)
 
  Output: polyphen_fpr.txt #(Final PolyPhen-2 result file)
 
  Output: newest_Mutation_merged_withoutPolyphen_scores.txt #(Final TUSON input file with NA PolyPhen-2 scores)

13. Entropy of missense calculation
```
cd step13_correct_entropy_out_file
tail -n +2 ../step12_TUSON_input_prepare/TUSON_input_new_revised.txt | sort -t$'\t' -k 1 > TUSON_input_new_sorted.txt
```
  Input: TUSON_input_new_revised.txt #(from Step 12)
 
  Output: TUSON_input_new_sorted.txt #(sorted TUSON input file)
```
python entropy_missense.py TUSON_input_new_sorted.txt entropy_output.txt
```
  Input: TUSON_input_new_sorted.txt #(from this step)
 
  Output: entropy_output.txt #(Entropy of missense result)
```
awk -F'\t' '{print $1}' entropy_output1.txt > gene_set_sorted.txt
``` 
  Input: entropy_output1.txt #(from this step)
 
  Output: gene_set_sorted.txt #(used by TUSON)

```
R --vanilla --slave --args entropy_output1.txt ../../Gene_annotation/Step2_gene_merge/Gene_set_new.txt entropy_pvalue_output_file1 < p_value_entropy.R
```
  Input: entropy_output.txt #(from last step)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)
 
  Output: entropy_pvalue_output_file1 #(Generation of missense entropy P-value file with new format as TUSON input)
```
javac step12_correct_entropy_out_file.java #(correct format of results of Entropy of missense calculation)
java step12_correct_entropy_out_file
cd ..
```
  Input: entropy_pvalue_output_file1 #(from last step)
 
  Output: entropy_pvalue_output_corrected.txt #(Obtain final entropy information)

14. TUSON (Time-consuming)
```
cd step14_TUSON
R --vanilla --slave --args ../step13_correct_entropy_out_file/TUSON_input_new_sorted.txt ../step13_correct_entropy_out_file/entropy_pvalue_output_corrected.txt ../../Gene_annotation/Step2_gene_merge/Gene_set_new.txt STOP.codon.Norm.Factor.txt TUSON_output.txt < TUSON_TSG_AND_OG.R
cd ..
```
  Input: TUSON_input_new_sorted.txt #(from last step)

  Input: entropy_pvalue_output_file_corrected.txt #(from last step)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)

  Input: STOP.codon.Norm.Factor.txt #(provided by TUSON)
 
  Output: TUSON_output_newest.txt #(Mutation feature quantification defined by TUSON and gene ranking results)

15. MGAentropy and ExonCons feature generation
```
cd step15_SNVBOX_processing
unzip snvbox.zip
cd ..
javac step15_SNVbox_features.java
java step15_SNVbox_features
cd ..
```
  Input: SNVBOX_results folder generated from CRAVAT website (http://hg19.cravat.us/CRAVAT/ with SnvGet options)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)

  Input: Ensembl_TransID_v87_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)

  Input: Ensembl_TransID_v92_GeneSymbol.txt #(Annotation file: Gene Symbol to Ensembl IDs for some outdated genes)

  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)

  Input: CosmicTranscripts.tsv #(from step 4)
  
  Input: Cosmic2Genename.txt #(Mapping file from COSMIC Gene ID to gene symbol)

  Internal: MGAentropy_stats.txt #(Internal result)
  
  Internal: Exoncons_stats.txt #(Internal result)
 
  Output: SNVBox_features.txt

16. Additional mutation features
```
cd step16_Additional_mutation_features
unzip gencode.v25lift37.annotation.gtf.zip
javac step16_Additional_mutation_features.java
java step16_Additional_mutation_features
cd ..
```
  Input: TUSON_input_new_sorted.txt #(from step 12)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)

  Input: CosmicMutantExport.tsv #(Cosmic Mutation file downloaded from COSMIC)

  Input: Ensembl_CDS_size.txt #(CDS information from Ensembl)

  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)

  Input: gencode.v25lift37.annotation.gtf #(Gencode v25 hg19 GTF annotation)

  Input: compiled_genelist.txt #(From step 2 from Gene_annotation folder)

  Input: compiled_genelist_CDS_length.txt #(generated from this step)
 
  Output: additional_mutation_features.txt #(Several additional features)
  
17. gnomAD data set

```
cd step17_gnomad_features
javac step17_gnomad_features.java 
java step17_gnomad_features
cd ..
```
  Input: gnomad.v2.1.1.lof_metrics.by_gene.txt #(The gene-centric gnomAD dataset, https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz)

  Input: Gene_set_new.txt #(Gene annotation information used by this study)

  Input: ensembl_synonymous_to_symbol.txt #(Gene alias conversion to HGNC gene names)
 
  Output: gnomad_LOF_processed.txt #(features in the gnomAD dataset)