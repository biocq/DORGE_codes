import java.io.*;
import java.math.BigDecimal;
import java.util.*;

public class github_feature_combination {
    public static void main(String[] args) throws Exception {
        String machine = "mac";
        combine_different_features(machine);
    }
   
	public static void combine_different_features(String machine) throws Exception {

        String path = "";
        String Line = "";
        String[] temp;

        if (machine.equals("mac")) {
            path = "";
        } else {
            path = "";
        }

        Map<String, String> Silent_KB_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> N_Missense = new LinkedHashMap<String, String>();
        Map<String, String> N_LOF = new LinkedHashMap<String, String>();
        Map<String, String> N_Splice = new LinkedHashMap<String, String>();

        Map<String, String> Missense_KB_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> LOF_KB_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Missense_entropy = new LinkedHashMap<String, String>();
        Map<String, String> Accurate_Missense_entropy = new LinkedHashMap<String, String>();
        Map<String, String> LOF_TO_Silent_Ratio = new LinkedHashMap<String, String>();

        Map<String, String> Splice_TO_Silent_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> nonSilent_TO_Silent_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Missense_TO_Silent_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Missense_Damaging_TO_Missense_Benign_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> LOF_TO_Benign_Ratio = new LinkedHashMap<String, String>();

        Map<String, String> Splice_TO_Benign_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Missense_TO_Benign_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Missense_Damaging_TO_Benign_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Polyphen2 = new LinkedHashMap<String, String>();

        Map<String, String> LOF_TO_Total_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Missense_TO_Total_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> Splice_TO_Total_Ratio = new LinkedHashMap<String, String>();
        Map<String, String> LOF_TO_Missense_Ratio = new LinkedHashMap<String, String>();

        Map<String, String> Silent_fraction = new LinkedHashMap<String, String>();
        Map<String, String> nonsense_fraction = new LinkedHashMap<String, String>();
        Map<String, String> Splice_site_fraction = new LinkedHashMap<String, String>();
        Map<String, String> Missense_fraction = new LinkedHashMap<String, String>();
        Map<String, String> Recurrent_missense_fraction = new LinkedHashMap<String, String>();
        Map<String, String> Frameshift_indel_fraction = new LinkedHashMap<String, String>();
        Map<String, String> Inframe_indel_fraction = new LinkedHashMap<String, String>();
        Map<String, String> Lost_start_and_stop_fraction = new LinkedHashMap<String, String>();
        Map<String, String> inactivating_mutations_fraction = new LinkedHashMap<String, String>();
        Map<String, String> CDS_length = new LinkedHashMap<String, String>();

        Map<String, String> log_gene_length = new LinkedHashMap<String, String>();

        BufferedReader brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/replication_time/compiled_gene_pos_hg19.bed"));//gene position
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            int length1=Integer.parseInt(temp[2])-Integer.parseInt(temp[1]);
            log_gene_length.put(temp[3].toUpperCase(),""+div(Math.log(length1),Math.log(2),2));
        }
        brc1.close();

        brc1 = new BufferedReader(new FileReader(
                path + "Mutation_processing/step14_TUSON/TUSON_output.txt"));//output_of_TUSON_TSG_AND_OG_additional_features.R
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            if(temp[29].equals("NA")){
                temp[29]="0";
            }
            Silent_KB_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[29]),1,2));
            int temp1=0;

            temp1=Integer.parseInt(temp[9]);
            temp1=temp1+1;//Pseudo count

            N_Missense.put(temp[0].toUpperCase(), div(Math.log(temp1),Math.log(2),2)+"");
            temp1=Integer.parseInt(temp[4]);
            temp1=temp1+1;//Pseudo count
            N_LOF.put(temp[0].toUpperCase(), div(Math.log(temp1),Math.log(2),2)+"");
            temp1=Integer.parseInt(temp[5]);
            temp1=temp1+1;//Pseudo count
            N_Splice.put(temp[0].toUpperCase(), div(Math.log(temp1),Math.log(2),2)+"");
            if(temp[30].equals("NA")){
                temp[30]="0";
            }
            Missense_KB_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[30]),1,2));
            if(temp[31].equals("NA")){
                temp[31]="0";
            }
            LOF_KB_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[31]),1,2));
            if(temp[42].equals("NA")){
                temp[42]="0";
            }
            if(temp[44].equals("NA")){
                temp[44]="0";
            }
            Missense_entropy.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[44]),1,2));
            Accurate_Missense_entropy.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[44]),1,6));
            LOF_TO_Silent_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[16]),1,2));


            double silent_mut=Double.parseDouble(temp[1]);
            silent_mut=silent_mut+1;//Pseudo count
            nonSilent_TO_Silent_Ratio.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[15])-silent_mut,silent_mut,2));
            Splice_TO_Silent_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[17]),1,2));
            Missense_TO_Silent_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[18]),1,2));
            Missense_Damaging_TO_Missense_Benign_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[19]),1,2));
            LOF_TO_Benign_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[20]),1,2));

            Splice_TO_Benign_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[21]),1,2));
            Missense_TO_Benign_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[22]),1,2));
            Missense_Damaging_TO_Benign_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[23]),1,2));
            if(temp[10].equals("NA")){
                temp[10]="0";
            }
            Polyphen2.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[10]),1,2));

            LOF_TO_Total_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[24]),1,2));
            Missense_TO_Total_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[25]),1,2));
            Splice_TO_Total_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[26]),1,2));
            LOF_TO_Missense_Ratio.put(temp[0].toUpperCase(), ""+div(Double.parseDouble(temp[27]),1,2));

        }
        brc1.close();

        brc1 = new BufferedReader(new FileReader(
                path + "Mutation_processing/step16_Additional_mutation_features/additional_mutation_features.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");

            if(temp[1].equals("NA")){
                temp[1]="0";
            }
            Silent_fraction.put(temp[0],""+div(Double.parseDouble(temp[1]),1,2));
            if(temp[2].equals("NA")){
                temp[2]="0";
            }
            nonsense_fraction.put(temp[0],""+div(Double.parseDouble(temp[2]),1,2));
            if(temp[3].equals("NA")){
                temp[3]="0";
            }
            Splice_site_fraction.put(temp[0],""+div(Double.parseDouble(temp[3]),1,2));
            if(temp[4].equals("NA")){
                temp[4]="0";
            }
            Missense_fraction.put(temp[0],""+div(Double.parseDouble(temp[4]),1,2));
            if(temp[5].equals("NA")){
                temp[5]="0";
            }
            Recurrent_missense_fraction.put(temp[0],""+div(Double.parseDouble(temp[5]),1,2));
            if(temp[6].equals("NA")){
                temp[6]="0";
            }
            Frameshift_indel_fraction.put(temp[0],""+div(Double.parseDouble(temp[6]),1,2));
            if(temp[7].equals("NA")){
                temp[7]="0";
            }
            Inframe_indel_fraction.put(temp[0],""+div(Double.parseDouble(temp[7]),1,2));
            if(temp[8].equals("NA")){
                temp[8]="0";
            }
            Lost_start_and_stop_fraction.put(temp[0],""+div(Double.parseDouble(temp[8]),1,2));
            if(temp[9].equals("NA")){
                temp[9]="0";
            }
            inactivating_mutations_fraction.put(temp[0],""+div(Double.parseDouble(temp[9]),1,2));
            if(temp[10].equals("NA")){
                temp[10]="0";
            }
            CDS_length.put(temp[0],""+div(Double.parseDouble(temp[10]),1,2));
        }
        brc1.close();

        Map<String, String> length_H3K4me3 = new LinkedHashMap<String, String>();
        Map<String, String> broad_H3K4me3_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K4me3_height = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K4me3_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            length_H3K4me3.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K4me3_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K4me3_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K4me2_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K4me2_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K4me2_width = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K4me2_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K4me2_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K4me2_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K4me2_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K4me1_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K4me1_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K4me1_width = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K4me1_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K4me1_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K4me1_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K4me1_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K9ac_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K9ac_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K9ac_width = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K9ac_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K9ac_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K9ac_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K9ac_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K9me2_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K9me2_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K9me2_width = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K9me2_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K9me2_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K9me2_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K9me2_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K9me3_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K9me3_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K9me3_width = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K9me3_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K9me3_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K9me3_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K9me3_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K27ac_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K27ac_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K27ac_width = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K27ac_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K27ac_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K27ac_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K27ac_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K27me3_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K27me3_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K27me3_width = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K27me3_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K27me3_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K27me3_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K27me3_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K36me3_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K36me3_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K36me3_width = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K36me3_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K36me3_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K36me3_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K36me3_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H3K79me2_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H3K79me2_height = new LinkedHashMap<String, String>();
        Map<String, String> H3K79me2_width = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H3K79me2_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H3K79me2_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H3K79me2_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H3K79me2_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> broad_H4K20me1_percentage = new LinkedHashMap<String, String>();
        Map<String, String> H4K20me1_height = new LinkedHashMap<String, String>();
        Map<String, String> H4K20me1_width = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/histone_modification/ENCODE_H4K20me1_peak_info.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            H4K20me1_width.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            broad_H4K20me1_percentage.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
            H4K20me1_height.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[3]),1,2));
        }
        brc1.close();

        Map<String, String> CNA_amplification = new LinkedHashMap<String, String>();
        Map<String, String> CNA_deletion = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Genomics_processing/CNA/CNA_freq.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            CNA_amplification.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));

            CNA_deletion.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[2]),1,2));
        }
        brc1.close();

        Map<String, String> CRISPR_screening = new LinkedHashMap<String, String>();
        Map<String, String> minus_CRISPR_screening = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Phenotype_processing/CRISPR_screen/CRISPR_screening_dependency_scores.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");

            if(temp[1].contains("?")){
                CRISPR_screening.put(temp[0].toUpperCase(),"-0.15");//imputation by mean
                minus_CRISPR_screening.put(temp[0].toUpperCase(),"0.15");//imputation by mean
            }else{
                CRISPR_screening.put(temp[0].toUpperCase(),div(Double.parseDouble(temp[1]),1,2)+"");
                minus_CRISPR_screening.put(temp[0].toUpperCase(),div(Double.parseDouble(temp[2]),1,2)+"");
            }

        }
        brc1.close();

        Map<String, String> Canyon_genebody_hypermethylation = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/genebody_canyon/Canyon_methylation.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            double ratio=0;
            if(temp[1].contains("?")){
                ratio=0;
            }else{
                if(Double.parseDouble(temp[1])>0){
                    ratio=div(Double.parseDouble(temp[2]),
                            Double.parseDouble(temp[1]),2);
                }else{
                    ratio=div(Double.parseDouble(temp[1]),0.1,2);
                }
                if(ratio<0.001){
                    ratio=1;
                }

            }
            Canyon_genebody_hypermethylation.put(temp[0].toUpperCase(),ratio+"");

        }
        brc1.close();

        Map<String, String> promoter_differential_methyl = new LinkedHashMap<String, String>();
        Map<String, String> promoter_differential_methyl_minus = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/methylation/promoter_differential_methyl_COSMIC.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            double diff=0;
            if(temp[1].contains("?")){
                diff=0;
            }else {
                diff=div(Double.parseDouble(temp[1])-Double.parseDouble(temp[2]),1,2);
            }
            promoter_differential_methyl.put(temp[0].toUpperCase(),diff+"");
            promoter_differential_methyl_minus.put(temp[0].toUpperCase(),(0-diff)+"");
        }
        brc1.close();

        Map<String, String> genebody_differential_methyl = new LinkedHashMap<String, String>();
        Map<String, String> genebody_differential_methyl_minus = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/methylation/genebody_differential_methyl_COSMIC.txt"));//

        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            double diff=0;
            if(temp[1].contains("?")){
                diff=0;
            }else {
                diff=div(Double.parseDouble(temp[1])-Double.parseDouble(temp[2]),1,2);
            }

            genebody_differential_methyl.put(temp[0].toUpperCase(),diff+"");
            genebody_differential_methyl_minus.put(temp[0].toUpperCase(),(0-diff)+"");

        }
        brc1.close();

        Map<String, String> dysregulation = new LinkedHashMap<String, String>();
        Map<String, String> minus_dysregulation = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Phenotype_processing/Gene_expression/gene_expression_quant.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            if(temp[1].equals("?")){
                temp[1]="0";
            }
            double exp = div(Double.parseDouble(temp[1]),1,2);
            if(exp > 10){
                exp= 10;
            }
            if(exp < -10){
                exp = -10;
            }
            dysregulation.put(temp[0].toUpperCase(),""+exp);
            minus_dysregulation.put(temp[0].toUpperCase(),""+(0-exp));

        }
        brc1.close();

        Map<String, String> BioGRID_log_degree = new LinkedHashMap<String, String>();
        Map<String, String> BioGRID_clossness = new LinkedHashMap<String, String>();
        Map<String, String> BioGRID_betweenness = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Evaluation_processing/BIOGRID/BioGRID_network_node_metrics_processed.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {

            temp = Line.split("\t");
            if(temp[3].equals("0")){
                temp[3]="1";
            }
            double logdegree =  div(Math.log(Integer.parseInt(temp[3])),Math.log(2),4);

            BioGRID_log_degree.put(temp[0],""+div(logdegree,1,2));

            BioGRID_clossness.put(temp[0],temp[2]);
            BioGRID_betweenness.put(temp[0],""+temp[1]);

        }
        brc1.close();


        Map<String, String> S50_replication_timing = new LinkedHashMap<String, String>();
        Map<String, String> minus_S50_replication_timing = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/replication_time/Replication_timing_percent.txt"));//

        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            if(temp[1].equals("NA")){
                temp[1]="0.46";
            }
            S50_replication_timing.put(temp[0].toUpperCase(),""+div(Double.parseDouble(temp[1]),1,2));
            minus_S50_replication_timing.put(temp[0].toUpperCase(),""+(1-div(Double.parseDouble(temp[1]),1,2)));
        }
        brc1.close();

        Map<String, String> MGAentropy = new LinkedHashMap<String, String>();
        Map<String, String> Exon_Cons = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Mutation_processing/step15_SNVBOX_processing/SNVBox_features.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {

            temp = Line.split("\t");
            MGAentropy.put(temp[0],""+div(Double.parseDouble(temp[1]),1,2));
            Exon_Cons.put(temp[0],""+div(Double.parseDouble(temp[3]),1,2));

        }
        brc1.close();

        Map<String, String> VEST = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Phenotype_processing/VEST/VEST_values.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {

            temp = Line.split("\t");
            if(temp[1].equals("NaN")){
                VEST.put(temp[0],"0");
            }else{
                double avg=div(Double.parseDouble(temp[1]),1,2);
                VEST.put(temp[0],""+avg);
            }
        }
        brc1.close();

        Map<String, String> SE = new LinkedHashMap<String, String>();
        brc1 = new BufferedReader(new FileReader(
                path + "Epigenetics_processing/Super_enhancer/super_enhancer_stats.txt"));//
        brc1.readLine();
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            double percentage=div(Integer.parseInt(temp[1]),99,2);
            SE.put(temp[0],""+percentage);
        }
        brc1.close();

        Map<String, String> pLI = new LinkedHashMap<String, String>();
        Map<String, String> pRec = new LinkedHashMap<String, String>();
        Map<String, String> pNull = new LinkedHashMap<String, String>();

        Map<String, String> syn_z = new LinkedHashMap<String, String>();
        Map<String, String> mis_z = new LinkedHashMap<String, String>();
        Map<String, String> lof_z = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Mutation_processing/step17_gnomad_features/gnomad_LOF_processed.txt"));//
        brc1.readLine();

        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            pLI.put(temp[0],""+div(Double.parseDouble(temp[1]),1,2));
            pRec.put(temp[0],""+div(Double.parseDouble(temp[2]),1,2));
            pNull.put(temp[0],""+div(Double.parseDouble(temp[3]),1,2));
            syn_z.put(temp[0],""+div(Double.parseDouble(temp[4]),1,2));
            mis_z.put(temp[0],""+div(Double.parseDouble(temp[5]),1,2));
            lof_z.put(temp[0],""+div(Double.parseDouble(temp[6]),1,2));
        }
        brc1.close();

        Map<String, String> familyMemberCount = new LinkedHashMap<String, String>();
        Map<String, String> ncRVIS = new LinkedHashMap<String, String>();
        Map<String, String> ncGERP = new LinkedHashMap<String, String>();
        Map<String, String> RVIS_percentile = new LinkedHashMap<String, String>();
        Map<String, String> dn_to_ds_ratio = new LinkedHashMap<String, String>();
        Map<String, String> GDI = new LinkedHashMap<String, String>();
        Map<String, String> gene_age = new LinkedHashMap<String, String>();

        brc1 = new BufferedReader(new FileReader(
                path + "Genomics_processing/NCBOOST_features/NCBoost_geneDB_processed.txt"));//
        brc1.readLine();

        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            familyMemberCount.put(temp[0],""+div(Double.parseDouble(temp[1]),1,2));
            ncRVIS.put(temp[0],""+div(Double.parseDouble(temp[2]),1,2));
            ncGERP.put(temp[0],""+div(Double.parseDouble(temp[3]),1,2));
            RVIS_percentile.put(temp[0],""+div(Double.parseDouble(temp[6]),1,2));
            dn_to_ds_ratio.put(temp[0],""+div(Double.parseDouble(temp[4]),1,2));
            GDI.put(temp[0],""+div(Double.parseDouble(temp[7]),1,2));
            gene_age.put(temp[0],""+div(Double.parseDouble(temp[5]),1,2));
        }
        brc1.close();

        BufferedReader brc3 = new BufferedReader(new FileReader(path + "Gene_annotation/Step2_gene_merge/Gene_set_new.txt"));
        brc3.readLine();
        int index = 1;
        Map<String,String> gene_name_mapping = new LinkedHashMap<String,String>();
        while ((Line = brc3.readLine()) != null) {
            temp = Line.split("\t");
            gene_name_mapping.put(index+"",temp[0].toUpperCase());
            index++;
        }
        brc3.close();

        BufferedWriter bw = new BufferedWriter(new FileWriter(
                path+"Features/All_features_TSG.csv"));
        BufferedWriter bw1 = new BufferedWriter(new FileWriter(
                path+"Features/All_features_OG.csv"));
        BufferedWriter bw2 = new BufferedWriter(new FileWriter(
                path+"Features/All_features_NG.csv"));
        BufferedWriter bw3 = new BufferedWriter(new FileWriter(
                path+"Features/All_features.csv"));
        BufferedWriter bw4 = new BufferedWriter(new FileWriter(
                path+"../DORGE_codes/Figure_3/data/Missense_entropy_accurate_version.txt"));
        bw4.write("Gene\tMissense_entropy\n");

        brc1 = new BufferedReader(new FileReader(
                path + "Gene_annotation/Step2_gene_merge/Gene_set_new.txt"));
        brc1.readLine();
        String sep=",";
        String title="Gene"
                +sep+"Silent_mutations_kb"
                +sep+"log_Total_N_missense_mutations"
                +sep+"log_Total_N_LoF_mutations"
                +sep+"log_Total_N_of_splicing_mutations"
                +sep+"Missense_mutations_kb"
                +sep+"LoF_mutations_kb"
                +sep+"Missense_entropy"
                +sep+"LOF_silent_ratio"
                +sep+"Splice_silent_ratio"
                +sep+"Missense_silent_ratio"
                +sep+"HiFI_missense_LoFI_missense_ratio"
                +sep+"LOF_benign_ratio"
                +sep+"Splice_benign_ratio"
                +sep+"Missense_benign_ratio"
                +sep+"Missense_damaging_benign_ratio"
                +sep+"PolyPhen_2_score"
                +sep+"LOF_total_ratio"
                +sep+"Missense_total_ratio"
                +sep+"Splice_total_ratio"
                +sep+"LOF_missense_ratio"
                +sep+"Silent_fraction"
                +sep+"Nonsense_fraction"
                +sep+"Missense_fraction"
                +sep+"Recurrent_missense_fraction"
                +sep+"Frameshift_indel_fraction"
                +sep+"Inframe_indel_fraction"
                +sep+"Lost_start_and_stop_fraction"
                +sep+"Inactivating_fraction"
                +sep+"NonSilent_silent_ratio"
                +sep+"log_gene_length"
                +sep+"log_CDS_length"
                +sep+"CNA_deletion_percentage"
                +sep+"CNA_amplification_percentage"
                +sep+"Exon_conservation_phastCons_score"
                +sep+"Missense_MGAentropy"
                +sep+"Early_replication_timing"
                +sep+"Late_replication_timing"
                +sep+"VEST_score"
                +sep+"Gene_expression_Z_score"
                +sep+"Gene_expression_minus_of_Z_score"
                +sep+"Increase_of_cell_proliferation_by_CRISPR_Knock_down"
                +sep+"Decrease_of_cell_proliferation_by_CRISPR_Knock_down"
                +sep+"Super_enhancer_percentage"
                +sep+"BioGRID_betweenness"
                +sep+"BioGRID_clossness"
                +sep+"Log_BioGRID_degree"
                +sep+"Promoter_hypermethylation_in_cancer"
                +sep+"Promoter_hypomethylation_in_cancer"
                +sep+"Gene_body_hypermethylation_in_cancer"
                +sep+"Gene_body_hypomethylation_in_cancer"
                +sep+"Gene_body_canyon_hypermethylation_in_cancer"
                +sep+"pLI_score"
                +sep+"pRec_score"
                +sep+"pNull_score"
                +sep+"Synonymous_o_e_constraint"
                +sep+"Missense_o_e_constraint"
                +sep+"LoF_o_e_constraint"
                +sep+"Primate_dN_dS_ratio"
                +sep+"Gene_damage_index"
                +sep+"RVIS_percentile"
                +sep+"ncRVIS_score"
                +sep+"ncGERP_score"
                +sep+"Gene_age"
                +sep+"Family_member_count"
                +sep+"H3K4me3_peak_length"
                +sep+"Percentage_of_broad_H3K4me3_peaks"
                +sep+"Height_of_H3K4me3_peaks"
                +sep+"H3K4me2_peak_length"
                +sep+"Percentage_of_broad_H3K4me2_peaks"
                +sep+"Height_of_H3K4me2_peaks"
                +sep+"H3K4me1_peak_length"
                +sep+"Percentage_of_broad_H3K4me1_peaks"
                +sep+"Height_of_H3K4me1_peaks"
                +sep+"H3K36me3_peak_length"
                +sep+"Percentage_of_broad_H3K36me3_peaks"
                +sep+"Height_of_H3K36me3_peaks"
                +sep+"H3K27ac_peak_length"
                +sep+"Percentage_of_broad_H3K27ac_peaks"
                +sep+"Height_of_H3K27ac_peaks"
                +sep+"H3K27me3_peak_length"
                +sep+"Percentage_of_broad_H3K27me3_peaks"
                +sep+"Height_of_H3K27me3_peaks"
                +sep+"H3K9me3_peak_length"
                +sep+"Percentage_of_broad_H3K9me3_peaks"
                +sep+"Height_of_H3K9me3_peaks"
                +sep+"H3K9ac_peak_length"
                +sep+"Percentage_of_broad_H3K9ac_peaks"
                +sep+"Height_of_H3K9ac_peaks"
                +sep+"H3K9me2_peak_length"
                +sep+"Percentage_of_broad_H3K9me2_peaks"
                +sep+"Height_of_H3K9me2_peaks"
                +sep+"H3K79me2_peak_length"
                +sep+"Percentage_of_broad_H3K79me2_peaks"
                +sep+"Height_of_H3K79me2_peaks"
                +sep+"H4K20me1_peak_length"
                +sep+"Percentage_of_broad_H4K20me1_peaks"
                +sep+"Height_of_H4K20me1_peaks";

        bw.write(title+"\n");
        bw1.write(title+"\n");
        bw2.write(title+"\n");
        bw3.write(title+"\n");
        String nullvalue="0";
        while ((Line = brc1.readLine()) != null) {
            temp = Line.split("\t");
            String line="";
            line=temp[0];

						if(Accurate_Missense_entropy.containsKey(temp[0].toUpperCase())){//7
                bw4.write(line+"\t"+Accurate_Missense_entropy.get(temp[0].toUpperCase())+"\n");
            }else{
                bw4.write(line+"\t0\n");
            }
            if(Silent_KB_Ratio.containsKey(temp[0].toUpperCase())){//1
                line=line+sep+Silent_KB_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(N_Missense.containsKey(temp[0].toUpperCase())){//2
                line=line+sep+N_Missense.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(N_LOF.containsKey(temp[0].toUpperCase())){//3
                line=line+sep+N_LOF.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(N_Splice.containsKey(temp[0].toUpperCase())){//4
                line=line+sep+N_Splice.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(Missense_KB_Ratio.containsKey(temp[0].toUpperCase())){//5
                line=line+sep+Missense_KB_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(LOF_KB_Ratio.containsKey(temp[0].toUpperCase())){//6
                line=line+sep+LOF_KB_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Missense_entropy.containsKey(temp[0].toUpperCase())){//7
                line=line+sep+Missense_entropy.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(LOF_TO_Silent_Ratio.containsKey(temp[0].toUpperCase())){//8
                line=line+sep+LOF_TO_Silent_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Splice_TO_Silent_Ratio.containsKey(temp[0].toUpperCase())){//9
                line=line+sep+Splice_TO_Silent_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(Missense_TO_Silent_Ratio.containsKey(temp[0].toUpperCase())){//10
                line=line+sep+Missense_TO_Silent_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Missense_Damaging_TO_Missense_Benign_Ratio.containsKey(temp[0].toUpperCase())){//11
                line=line+sep+Missense_Damaging_TO_Missense_Benign_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(LOF_TO_Benign_Ratio.containsKey(temp[0].toUpperCase())){//12
                line=line+sep+LOF_TO_Benign_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Splice_TO_Benign_Ratio.containsKey(temp[0].toUpperCase())){//13
                line=line+sep+Splice_TO_Benign_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Missense_TO_Benign_Ratio.containsKey(temp[0].toUpperCase())){//14
                line=line+sep+Missense_TO_Benign_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Missense_Damaging_TO_Benign_Ratio.containsKey(temp[0].toUpperCase())){//15
                line=line+sep+Missense_Damaging_TO_Benign_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Polyphen2.containsKey(temp[0].toUpperCase())){//16
                line=line+sep+Polyphen2.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(LOF_TO_Total_Ratio.containsKey(temp[0].toUpperCase())){//17
                line=line+sep+LOF_TO_Total_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Missense_TO_Total_Ratio.containsKey(temp[0].toUpperCase())){//18
                line=line+sep+Missense_TO_Total_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Splice_TO_Total_Ratio.containsKey(temp[0].toUpperCase())){//19
                line=line+sep+Splice_TO_Total_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(LOF_TO_Missense_Ratio.containsKey(temp[0].toUpperCase())){//20
                line=line+sep+LOF_TO_Missense_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(Silent_fraction.containsKey(temp[0].toUpperCase())){//21
                line=line+sep+Silent_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(nonsense_fraction.containsKey(temp[0].toUpperCase())){//22
                line=line+sep+nonsense_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            //if(Splice_site_fraction.containsKey(temp[0].toUpperCase())){
            //    line=line+sep+Splice_site_fraction.get(temp[0].toUpperCase());
            //}else{
            //    line=line+sep+nullvalue;
            //}
            if(Missense_fraction.containsKey(temp[0].toUpperCase())){//23
                line=line+sep+Missense_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Recurrent_missense_fraction.containsKey(temp[0].toUpperCase())){//24
                line=line+sep+Recurrent_missense_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Frameshift_indel_fraction.containsKey(temp[0].toUpperCase())){//25
                line=line+sep+Frameshift_indel_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Inframe_indel_fraction.containsKey(temp[0].toUpperCase())){//26
                line=line+sep+Inframe_indel_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Lost_start_and_stop_fraction.containsKey(temp[0].toUpperCase())){//27
                line=line+sep+Lost_start_and_stop_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(inactivating_mutations_fraction.containsKey(temp[0].toUpperCase())){//28
                line=line+sep+inactivating_mutations_fraction.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(nonSilent_TO_Silent_Ratio.containsKey(temp[0].toUpperCase())){//29
                line=line+sep+nonSilent_TO_Silent_Ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(log_gene_length.containsKey(temp[0].toUpperCase())){//30
                line=line+sep+log_gene_length.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(CDS_length.containsKey(temp[0].toUpperCase())){//31
                line=line+sep+CDS_length.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(CNA_deletion.containsKey(temp[0].toUpperCase())){//32
                line=line+sep+CNA_deletion.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(CNA_amplification.containsKey(temp[0].toUpperCase())){//33
                line=line+sep+CNA_amplification.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Exon_Cons.containsKey(temp[0].toUpperCase())){//34
                line=line+sep+Exon_Cons.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(MGAentropy.containsKey(temp[0].toUpperCase())){//35
                line=line+sep+MGAentropy.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(S50_replication_timing.containsKey(temp[0].toUpperCase())){//36
                line=line+sep+S50_replication_timing.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(S50_replication_timing.containsKey(temp[0].toUpperCase())){//37
                line=line+sep+(1-Double.parseDouble(S50_replication_timing.get(temp[0].toUpperCase())));
            }else{
                line=line+sep+nullvalue;
            }
            if(VEST.containsKey(temp[0].toUpperCase())){//38
                line=line+sep+VEST.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(dysregulation.containsKey(temp[0].toUpperCase())){//39
                line=line+sep+dysregulation.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(minus_dysregulation.containsKey(temp[0].toUpperCase())){//40
                line=line+sep+minus_dysregulation.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(CRISPR_screening.containsKey(temp[0].toUpperCase())){//41
                line=line+sep+CRISPR_screening.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(minus_CRISPR_screening.containsKey(temp[0].toUpperCase())){//42
                line=line+sep+minus_CRISPR_screening.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(SE.containsKey(temp[0].toUpperCase())){//43
                line=line+sep+SE.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(BioGRID_betweenness.containsKey(temp[0].toUpperCase())){//44
                line=line+sep+BioGRID_betweenness.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(BioGRID_clossness.containsKey(temp[0].toUpperCase())){//45
                line=line+sep+BioGRID_clossness.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(BioGRID_log_degree.containsKey(temp[0].toUpperCase())){//46
                line=line+sep+BioGRID_log_degree.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(promoter_differential_methyl.containsKey(temp[0].toUpperCase())){//47
                line=line+sep+promoter_differential_methyl.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(promoter_differential_methyl_minus.containsKey(temp[0].toUpperCase())){//48
                line=line+sep+promoter_differential_methyl_minus.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(genebody_differential_methyl.containsKey(temp[0].toUpperCase())){//49
                line=line+sep+genebody_differential_methyl.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(genebody_differential_methyl_minus.containsKey(temp[0].toUpperCase())){//50
                line=line+sep+genebody_differential_methyl_minus.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(Canyon_genebody_hypermethylation.containsKey(temp[0].toUpperCase())){//51
                line=line+sep+Canyon_genebody_hypermethylation.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(pLI.containsKey(temp[0].toUpperCase())){//52
                line=line+sep+pLI.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(pRec.containsKey(temp[0].toUpperCase())){//53
                line=line+sep+pRec.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(pNull.containsKey(temp[0].toUpperCase())){//54
                line=line+sep+pNull.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(syn_z.containsKey(temp[0].toUpperCase())){//55
                line=line+sep+syn_z.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(mis_z.containsKey(temp[0].toUpperCase())){//56
                line=line+sep+mis_z.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(lof_z.containsKey(temp[0].toUpperCase())){//57
                line=line+sep+lof_z.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(dn_to_ds_ratio.containsKey(temp[0].toUpperCase())){//58
                line=line+sep+dn_to_ds_ratio.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(GDI.containsKey(temp[0].toUpperCase())){//59
                line=line+sep+GDI.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(RVIS_percentile.containsKey(temp[0].toUpperCase())){//60
                line=line+sep+RVIS_percentile.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(ncRVIS.containsKey(temp[0].toUpperCase())){//61
                line=line+sep+ncRVIS.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(ncGERP.containsKey(temp[0].toUpperCase())){//62
                line=line+sep+ncGERP.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(gene_age.containsKey(temp[0].toUpperCase())){//63
                line=line+sep+gene_age.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(familyMemberCount.containsKey(temp[0].toUpperCase())){//64
                line=line+sep+familyMemberCount.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(length_H3K4me3.containsKey(temp[0].toUpperCase())){//65
                line=line+sep+length_H3K4me3.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K4me3_percentage.containsKey(temp[0].toUpperCase())){//66
                line=line+sep+broad_H3K4me3_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K4me3_height.containsKey(temp[0].toUpperCase())){//67
                line=line+sep+H3K4me3_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(H3K4me2_width.containsKey(temp[0].toUpperCase())){//68
                line=line+sep+H3K4me2_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K4me2_percentage.containsKey(temp[0].toUpperCase())){//69
                line=line+sep+broad_H3K4me2_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K4me2_height.containsKey(temp[0].toUpperCase())){//70
                line=line+sep+H3K4me2_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(H3K4me1_width.containsKey(temp[0].toUpperCase())){//71
                line=line+sep+H3K4me1_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K4me1_percentage.containsKey(temp[0].toUpperCase())){//72
                line=line+sep+broad_H3K4me1_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K4me1_height.containsKey(temp[0].toUpperCase())){//73
                line=line+sep+H3K4me1_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K36me3_width.containsKey(temp[0].toUpperCase())){//74
                line=line+sep+H3K36me3_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K36me3_percentage.containsKey(temp[0].toUpperCase())){//75
                line=line+sep+broad_H3K36me3_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K36me3_height.containsKey(temp[0].toUpperCase())){//76
                line=line+sep+H3K36me3_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K27ac_width.containsKey(temp[0].toUpperCase())){//77
                line=line+sep+H3K27ac_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K27ac_percentage.containsKey(temp[0].toUpperCase())){//78
                line=line+sep+broad_H3K27ac_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K27ac_height.containsKey(temp[0].toUpperCase())){//79
                line=line+sep+H3K27ac_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }

            if(H3K27me3_width.containsKey(temp[0].toUpperCase())){//80
                line=line+sep+H3K27me3_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K27me3_percentage.containsKey(temp[0].toUpperCase())){//81
                line=line+sep+broad_H3K27me3_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K27me3_height.containsKey(temp[0].toUpperCase())){//82
                line=line+sep+H3K27me3_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K9me3_width.containsKey(temp[0].toUpperCase())){//83
                line=line+sep+H3K9me3_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K9me3_percentage.containsKey(temp[0].toUpperCase())){//84
                line=line+sep+broad_H3K9me3_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K9me3_height.containsKey(temp[0].toUpperCase())){//85
                line=line+sep+H3K9me3_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K9ac_width.containsKey(temp[0].toUpperCase())){//86
                line=line+sep+H3K9ac_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K9ac_percentage.containsKey(temp[0].toUpperCase())){//87
                line=line+sep+broad_H3K9ac_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K9ac_height.containsKey(temp[0].toUpperCase())){//88
                line=line+sep+H3K9ac_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K9me2_width.containsKey(temp[0].toUpperCase())){//89
                line=line+sep+H3K9me2_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K9me2_percentage.containsKey(temp[0].toUpperCase())){//90
                line=line+sep+broad_H3K9me2_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K9me2_height.containsKey(temp[0].toUpperCase())){//91
                line=line+sep+H3K9me2_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K79me2_width.containsKey(temp[0].toUpperCase())){//92
                line=line+sep+H3K79me2_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H3K79me2_percentage.containsKey(temp[0].toUpperCase())){//93
                line=line+sep+broad_H3K79me2_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H3K79me2_height.containsKey(temp[0].toUpperCase())){//94
                line=line+sep+H3K79me2_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H4K20me1_width.containsKey(temp[0].toUpperCase())){//95
                line=line+sep+H4K20me1_width.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(broad_H4K20me1_percentage.containsKey(temp[0].toUpperCase())){//96
                line=line+sep+broad_H4K20me1_percentage.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }
            if(H4K20me1_height.containsKey(temp[0].toUpperCase())){//97
                line=line+sep+H4K20me1_height.get(temp[0].toUpperCase());
            }else{
                line=line+sep+nullvalue;
            }


            if (temp[4].equals("1")) {
                bw.write(line.replaceAll(",NA,",",0,").replaceAll(",,",",0.0,")+"\n");
            }
            if (temp[5].equals("1")) {
                bw1.write(line.replaceAll(",NA,",",0,").replaceAll(",,",",0.0,")+"\n");
            }
            if (temp[3].equals("1")) {
                bw2.write(line.replaceAll(",NA,",",0,").replaceAll(",,",",0.0,")+"\n");
            }
            bw3.write(line.replaceAll(",NA,",",0,").replaceAll(",,",",0.0,")+"\n");
        }

        brc1.close();
        bw.close();
        bw1.close();
        bw2.close();
        bw3.close();
        bw4.close();
    }


  	@SuppressWarnings("deprecation")
    public static double div(double v1, double v2, int scale) {
        if (scale < 0) {
            throw new IllegalArgumentException(
                    "The scale must be a positive integer or zero");
        }
        BigDecimal b1 = new BigDecimal(Double.toString(v1));
        BigDecimal b2 = new BigDecimal(Double.toString(v2));
        return b1.divide(b2, scale, BigDecimal.ROUND_HALF_UP).doubleValue();
    }

    public static ArrayList refreshFileList(String strPath,ArrayList filelist) {
        File dir = new File(strPath);

        File[] files = dir.listFiles();

        for (int i = 0; i < files.length; i++) {
            if (files[i].isDirectory()) {
                refreshFileList(files[i].getAbsolutePath(),filelist);
            } else {
                String strFileName = files[i].getAbsolutePath();
                filelist.add(strFileName);
            }
        }
        return filelist;
    }


}