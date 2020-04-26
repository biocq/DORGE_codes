import java.io.*;
import java.math.BigDecimal;
import java.util.*;


public class github_prediction_TSGOG_evaluation {

    public static void main(String[] args) throws Exception {
        String machine = "mac";

        
        different_approach_performance_TSGs_OGs(machine);//TSG/OG prediction evaluation - Data file S2
    }


    public static void different_approach_performance_TSGs_OGs(String machine) throws Exception {
        different_approach_performance_TSGs(machine,"TSG");
        different_approach_performance_TSGs(machine,"OG");
        different_approach_performance_TSGs_same_number_with_DORGE(machine,"TSG");
        different_approach_performance_TSGs_same_number_with_DORGE(machine,"OG");
    }


    public static void different_approach_performance_TSGs(String machine,String TSGorOG) throws Exception {
        String path = "";
        String Line = "";
        String[] temp;

        if (machine.equals("mac")) {
            path = "";
        } else {
            path = "";
        }
			int index=0;
        Map<String, String> genes = new LinkedHashMap<String, String>();

        BufferedWriter bw = new BufferedWriter(new FileWriter(
                path + ""+TSGorOG+"_gene_performance.txt"));
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
        
                
        String[] predicted=new String[genes.size()];
        String result="";
        brc = new BufferedReader(new FileReader(
                path+"../../DORGE_prediction.txt"));
        brc.readLine();
        while ((Line = brc.readLine()) != null) {
            temp = Line.split("\t");
            double TSG_prob=Double.parseDouble(temp[1]);
            double OG_prob=Double.parseDouble(temp[2]);
            
            if(TSGorOG.equals("TSG")){
            	if(TSG_prob>0.62485){
              		genes.put(temp[0],"");
            	}
		        }else{
		           if(OG_prob>0.7004394){
                	genes.put(temp[0],"");
            	}
		        }
        }
        brc.close();

        predicted=new String[genes.size()];
        index=0;
        Iterator<Map.Entry<String, String>> name1 = genes.entrySet().iterator();
        while (name1.hasNext()) {
            Map.Entry<String, String> entry = (Map.Entry<String, String>) name1
                    .next();
            String name = entry.getKey();
            predicted[index]=name;
            index++;
        }
        genes.clear();
        if(TSGorOG.equals("TSG")){
            result=performance_different_approaches_on_TSGs("DORGE",TSG,NG,all_genes,predicted);

        }else{
            result=performance_different_approaches_on_TSGs("DORGE",OG,NG,all_genes,predicted);
        }
        

        brc = new BufferedReader(new FileReader(
                path+"q0point1/2020_"+TSGorOG+".txt"));

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

        if(TSGorOG.equals("TSG")){
            result=result+performance_different_approaches_on_TSGs("20/20+",TSG,NG,all_genes,predicted);
        }else{
            result=result+performance_different_approaches_on_TSGs("20/20+",OG,NG,all_genes,predicted);
        }

        brc = new BufferedReader(new FileReader(
                path+"q0point1/TUSON_"+TSGorOG+"_reported_by_2020plus.txt"));
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
        if(TSGorOG.equals("TSG")){
            result=result+performance_different_approaches_on_TSGs("TUSON",TSG,NG,all_genes,predicted);

        }else{
            result=result+performance_different_approaches_on_TSGs("TUSON",OG,NG,all_genes,predicted);

        }
        
        
        brc = new BufferedReader(new FileReader(
                path+"q0point1/nonredudant_GUST_"+TSGorOG+".txt"));
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
        if(TSGorOG.equals("TSG")){
            result=result+performance_different_approaches_on_TSGs("GUST",TSG,NG,all_genes,predicted);

        }else{
            result=result+performance_different_approaches_on_TSGs("GUST",OG,NG,all_genes,predicted);

        }

        bw.write(result);

        bw.close();

    }

    public static void different_approach_performance_TSGs_same_number_with_DORGE(String machine,String TSGorOG) throws Exception {
        String path = "";
        String Line = "";
        String[] temp;

        if (machine.equals("mac")) {
            path = "";
        } else {
            path = "";
        }

        Map<String, String> genes = new LinkedHashMap<String, String>();

        BufferedWriter bw = new BufferedWriter(new FileWriter(
                path + ""+TSGorOG+"_gene_performance_same_gene_number_with_DORGE.txt"));
        Map<String, Boolean> TSG = new LinkedHashMap<String, Boolean>();
        Map<String, Boolean> OG = new LinkedHashMap<String, Boolean>();
        Map<String, Boolean> NG = new LinkedHashMap<String, Boolean>();
        Map<String, Boolean> all_genes = new LinkedHashMap<String, Boolean>();
        bw.write("Algorithm\t#\tSn\tSp\tPrecision\tAccuracy\n");
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
                path+"topgenes/2020_"+TSGorOG+".txt"));

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
        String result="";
        if(TSGorOG.equals("TSG")){
            result=performance_different_approaches_on_TSGs("2020",TSG,NG,all_genes,predicted);
        }else{
            result=performance_different_approaches_on_TSGs("2020",OG,NG,all_genes,predicted);
        }

        brc = new BufferedReader(new FileReader(
                path+"topgenes/TUSON_"+TSGorOG+".txt"));
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
        if(TSGorOG.equals("TSG")){
            result=result+performance_different_approaches_on_TSGs("TUSON",TSG,NG,all_genes,predicted);

        }else{
            result=result+performance_different_approaches_on_TSGs("TUSON",OG,NG,all_genes,predicted);

        }

        bw.write(result);

        bw.close();

    }


    static String performance_different_approaches_on_TSGs(String label,Map<String, Boolean> TSG
            ,Map<String, Boolean> NG,Map<String, Boolean> all_genes, String[] target)  {

        int tp = 0;
        int fp = 0;
        int tn = 0;
        int fn = 0;
        Map<String, Boolean> all_backup =new  LinkedHashMap<String, Boolean>();
        Iterator<Map.Entry<String, Boolean>> name1 = all_genes.entrySet().iterator();

        while (name1.hasNext()) {
            Map.Entry<String, Boolean> entry = (Map.Entry<String, Boolean>) name1
                    .next();
            String name = entry.getKey();
            all_backup.put(name,false);
        }
        for(int i=0;i<target.length;i++){
            if(TSG.containsKey(target[i].toUpperCase())){
                tp++;
                all_backup.remove(target[i].toUpperCase());
                //System.out.println(i+" "+target[i]);
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
        //all_genes=all_backup;
        for(int i=0;i<others.length;i++){
            if(TSG.containsKey(others[i].toUpperCase())){
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
