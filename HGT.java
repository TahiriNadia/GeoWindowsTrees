/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.File;

/**
 *
 * @author Nadia Tahiri
 */
public class HGT {

     public void Criterion(String fichier) throws IOException, InterruptedException{
       try {
           Runtime runtime = Runtime.getRuntime();
           Process process = runtime.exec("executable/criterion Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/"+fichier+" Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/output.txt Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/tmp Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/matrice.txt");
           ///Process process = runtime.exec("executable/criterion -inputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/"+fichier+" -outputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/output.txt");
		   InputStream is = process.getInputStream();
           InputStreamReader isr = new InputStreamReader(is);
           BufferedReader br = new BufferedReader(isr);
           String line;

           while ((line = br.readLine()) != null) {
             System.out.println(line);
           }

           process.waitFor();
       }catch (Exception e) {
       }
	   
    }	


     public boolean CriterionConsence(String fichier) throws IOException, InterruptedException{
       boolean error = false ;
       try {
           Runtime runtime = Runtime.getRuntime();
           Process process = runtime.exec("executable/criterion Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/"+fichier+" Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/outputConsence.txt Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/tmp Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/matrice.txt");
		   //Process process = runtime.exec("executable/criterion Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/"+fichier+" Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/outputConsence.txt tmp matrice.txt");
           //Process process = runtime.exec("executable/criterion -inputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/"+fichier+" -outputfile=Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/outputConsence.txt");
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
       
       File fichierConsense=new File("");
       StringBuilder st=new StringBuilder(); 
       Runtime runtime = Runtime.getRuntime();	
       File fi=new File("");
       //runtime.exec("cp "+fichierConsense.getAbsolutePath()+"/results.txt  "+fichierConsense.getAbsolutePath()+"/Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/outputConsence.txt");
       //runtime.exec("chmod +x "+fichierConsense.getAbsolutePath()+"/Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/outputConsence.txt");
       //System.out.println("cp "+fichierConsense.getAbsolutePath()+"/results.txt  "+fichierConsense.getAbsolutePath()+"/Result/Trees/"+fichier.substring(0,fichier.length()-4)+"/outputConsence.txt");

       if (error){
           return error;
       }else{
           /* System.out.print("Traitement terminé de Criterion Consence.\n"); */
           return error;
       }
    }
}
