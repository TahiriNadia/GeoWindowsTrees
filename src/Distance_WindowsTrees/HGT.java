/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 *
 * @author nadia
 */
public class HGT {

     public void Criterion(String fichier) throws IOException, InterruptedException{
       try {
           Runtime runtime = Runtime.getRuntime();
           Process process = runtime.exec("executable/criterion -inputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/"+fichier+" -outputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/output.txt");
           InputStream is = process.getInputStream();
           InputStreamReader isr = new InputStreamReader(is);
           BufferedReader br = new BufferedReader(isr);
           String line;

           while ((line = br.readLine()) != null) {
             System.out.println(line);
           }

           process.waitFor();
           System.out.print("Traitement terminé de Criterion.\n");
       }catch (Exception e) {
       }
    }


     public boolean CriterionConsence(String fichier) throws IOException, InterruptedException{
       boolean error = false ;
       try {
           Runtime runtime = Runtime.getRuntime();
           Process process = runtime.exec("executable/criterion -inputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/"+fichier+" -outputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/outputConsence.txt");
           InputStream is = process.getInputStream();
           InputStreamReader isr = new InputStreamReader(is);
           BufferedReader br = new BufferedReader(isr);
           
           String line;

           while ((line = br.readLine()) != null) {
               if (line.equalsIgnoreCase("Use valgrind or gdb to fix the problem")){
                   error = true;
               }
               System.out.println(line);
           }

           process.waitFor();
       }catch (Exception e) {
       }
       if (error){
           return error;
       }else{
           System.out.print("Traitement terminé de Criterion Consence.\n");
           return error;
       }
    }
}
