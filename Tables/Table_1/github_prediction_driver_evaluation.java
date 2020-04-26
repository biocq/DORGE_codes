import java.io.*;
import java.math.BigDecimal;
import java.util.*;

public class github_prediction_driver_evaluation {

    public static void main(String[] args) throws Exception {
        String machine = "mac";

        different_approach_performance(machine);//Driver gene prediction evaluation - Table 1
        
    }

    public static void different_approach_performance(String machine) throws Exception {
        String path = "";
        String Line = "";
        String[] temp;

        if (machine.equals("mac")) {
            path = "";
        } else {
            path = "";
        }

        Map<String, String> genes = new LinkedHashMap<String, String>();
        Map<String, Integer> count = new LinkedHashMap<String, Integer>();

        BufferedWriter bw = new BufferedWriter(new FileWriter(
                path + "Table1_driver_gene_performance.txt"));
        Map<String, Boolean> TSG = new LinkedHashMap<String, Boolean>();
        Map<String, Boolean> OG = new LinkedHashMap<String, Boolean>();
        Map<String, Boolean> NG = new LinkedHashMap<String, Boolean>();
        Map<String, Boolean> all_genes = new LinkedHashMap<String, Boolean>();
        bw.write("Method\t#\tSn\tSp\tPrecision\tAccuracy\n");
        BufferedReader brc = new BufferedReader(new FileReader(path + "../../Gene_set_new.txt"));
        brc.readLine();

        while ((Line = brc.readLine()) != null) {
            temp = Line.split("\t");
            
            all_genes.put(temp[0], false);
            if (temp[4].equals("1")) {
                TSG.put(temp[0], false);
            }
            if (temp[5].equals("1")) {
                OG.put(temp[0], false);
            }
            if (temp[3].equals("1")) {
                NG.put(temp[0], false);
            }
        }
        brc.close();


        brc = new BufferedReader(new FileReader(
                path+"q0point1/2020_driver.txt"));

        while ((Line = brc.readLine()) != null) {

            genes.put(Line,"");

        }
        brc.close();

        String[] predicted=new String[genes.size()];
        int index=0;
        Iterator<Map.Entry<String, String>> name1 = genes.entrySet().iterator();

        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        String result=performance_different_approaches_on_cancer_driver("20/20+",TSG,OG,NG,all_genes,predicted);

        brc = new BufferedReader(new FileReader(
                path+"../../DORGE_prediction.txt"));
        brc.readLine();
        while ((Line = brc.readLine()) != null) {
            temp = Line.split("\t");
            double TSG_prob=Double.parseDouble(temp[1]);
            double OG_prob=Double.parseDouble(temp[2]);
            if(TSG_prob>0.62485){
                genes.put(temp[0],"");
            }
            if(OG_prob>0.7004394){
                genes.put(temp[0],"");
            }

        }
        brc.close();
        

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("DORGE",TSG,OG,NG,all_genes,predicted);

        brc = new BufferedReader(new FileReader(
                path+"q0point1/ActiveDriver.txt"));
        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("ActiveDriver",TSG,OG,NG,all_genes,predicted);


  			brc = new BufferedReader(new FileReader(
                path+"q0point1/MutPanning.txt"));
        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("MutPanning",TSG,OG,NG,all_genes,predicted);
        
        brc = new BufferedReader(new FileReader(
                path+"q0point1/MuSIC.txt"));
        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("MuSIC",TSG,OG,NG,all_genes,predicted);

        brc = new BufferedReader(new FileReader(
                path+"q0point1/MutSigCV.txt"));
        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();
        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("MutSigCV",TSG,OG,NG,all_genes,predicted);

        brc = new BufferedReader(new FileReader(
                path+"q0point1/OncodriveCLUST.txt"));
        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("OncodriveCLUST",TSG,OG,NG,all_genes,predicted);

        brc = new BufferedReader(new FileReader(
                path+"q0point1/OncodriveFM.txt"));
        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();
        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("OncodriveFM",TSG,OG,NG,all_genes,predicted);

        brc = new BufferedReader(new FileReader(
                path+"q0point1/OncodriveFML.txt"));

        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("OncodriveFML",TSG,OG,NG,all_genes,predicted);

        brc = new BufferedReader(new FileReader(
                path+"q0point1/nonredudant_GUST_driver.txt"));
        while ((Line = brc.readLine()) != null) {
            genes.put(Line,"");
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("GUST",TSG,OG,NG,all_genes,predicted);


        brc = new BufferedReader(new FileReader(
                path+"q0point1/TUSON_driver_reported_by_2020plus.txt"));
        while ((Line = brc.readLine()) != null) {
            //temp = Line.split("\t");
            genes.put(Line,"");
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        result=result+performance_different_approaches_on_cancer_driver("TUSON",TSG,OG,NG,all_genes,predicted);
        bw.write(result);
        bw.close();

    }


    static String performance_different_approaches_on_cancer_driver(String label,Map<String, Boolean> TSG,Map<String, Boolean> OG
            ,Map<String, Boolean> NG,Map<String, Boolean> all_genes, String[] target)  {

        int tp = 0;
        int fp = 0;
        int tn = 0;
        int fn = 0;
        Map<String, Boolean> all_backup =new  LinkedHashMap<String, Boolean>() ;
        Iterator<Map.Entry<String, Boolean>> name1 = all_genes.entrySet().iterator();

        while (name1.hasNext()) {
            Map.Entry<String, Boolean> entry = (Map.Entry<String, Boolean>) name1
                    .next();
            String name = entry.getKey();
            all_backup.put(name,false);
        }
        for(int i=0;i<target.length;i++){
            if(TSG.containsKey(target[i].toUpperCase())||OG.containsKey(target[i].toUpperCase())){
                tp++;
                all_backup.remove(target[i].toUpperCase());
            }else if(NG.containsKey(target[i].toUpperCase())){
                fp++;
                all_backup.remove(target[i].toUpperCase());
            }else if(!all_backup.containsKey(target[i].toUpperCase())){
                all_backup.remove(target[i].toUpperCase());
            }
        }

        String[] others=new String[all_backup.size()];
        int index=0;
        name1 = all_backup.entrySet().iterator();

        while (name1.hasNext()) {
            Map.Entry<String, Boolean> entry = (Map.Entry<String, Boolean>) name1
                    .next();
            String name = entry.getKey();
            others[index]=name;
            index++;
        }

        for(int i=0;i<others.length;i++){
            if(TSG.containsKey(others[i].toUpperCase())||OG.containsKey(others[i].toUpperCase())){
                fn++;
            }else if(NG.containsKey(others[i].toUpperCase())){
                tn++;
            }
        }
        double sn=div(tp,tp+fn,3);//recall
        double sp=div(tn,fp+tn,3);
        double precision=div(tp,tp+fp,3);//ppv
        double accuracy=div(tp+tn,tp+fp+tn+fn,3);
        //System.out.println(label+"\t"+tp+"\t"+fp+"\t"+tn+"\t"+fn+"\n");
        return(label+"\t"+target.length+"\t"+sn+"\t"+sp+"\t"+precision+"\t"+accuracy+"\n");

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

}