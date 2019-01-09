/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Scanner;

/**
 *
 * @author Nadia Tahiri
 */
public class Phylip {


    public static String current=(new File("").getAbsolutePath());
    public  File repertoire=new File(current);


    /**
     * copie le fichier source dans le fichier resultat
     * retourne vrai si cela réussit
     */
    public static boolean Copy(boolean modeAjout, String nomDossier,String source, String dest, String pasWin){

        try{
            FileWriter destFile = new FileWriter(current+"/Result/Phylip/"+nomDossier+"/"+pasWin+"/"+dest,modeAjout);

            InputStream ips=new FileInputStream(source);
            InputStreamReader ipsr=new InputStreamReader(ips);
            BufferedReader br=new BufferedReader(ipsr);
            String ligne;

            while ((ligne=br.readLine())!=null){
                destFile.write(ligne);
                destFile.write("\n");
            }

            destFile.close();
            br.close();


        }catch (Exception e){
            return false; // Erreur
        }
       return true; // Résultat OK

    }


	public boolean createCommandFile(String filename, String commandeLine) {
        try {
            BufferedWriter command_file=new BufferedWriter(new FileWriter(new File(filename)));
            //File f=new File("");
	    //String command="#!/bin/sh\nexecutable/seqboot <"+f.getAbsolutePath()+"/executable/Inputs/inputSB >sortieSB\n";
	    String command=commandeLine;
           command_file.append(command);
           command_file.flush();
           command_file.close();
        } catch(Exception e) {
            e.printStackTrace();
		return false;
        }
        return true;
    }

       public void SeqBoot(String NameGene, String pasWin) throws IOException, InterruptedException{
	try {
	StringBuilder st=new StringBuilder(); 
        Runtime runtime = Runtime.getRuntime();	
	File fi=new File("");
	createCommandFile("executeSeqboot", "#!/bin/sh\nexecutable/seqboot <"+fi.getAbsolutePath()+"/executable/Inputs/inputSB >sortieSB\n");
        runtime.exec("chmod +x executeSeqboot");
	String[] command=new String[5];
        for (int i=0; i<command.length;i++) command[i]="";
        //Execute the Makefile_test.command.sh  
        File f=new File("");
        command[0]=f.getAbsolutePath()+"/executeSeqboot";
              
	       Process p = runtime.exec(command);
               InputStream stdoutput = p.getInputStream();
               InputStream stderr = p.getErrorStream();
               InputStreamReader isr = new InputStreamReader(stdoutput);
               InputStreamReader isr2 = new InputStreamReader(stderr);
               BufferedReader br = new BufferedReader(isr);
               BufferedReader br2 = new BufferedReader(isr2);
               String line = null;
               //--System.out
               while ( (line = br.readLine()) != null) {
                 st.append(line+"\n");
                  System.out.println(line);
               }
                   int exitVal = p.waitFor();                                                                
                  //--System.err
                   if (exitVal!=0) {
                       System.out.println("Seqboot error "+exitVal);
                      
                     while ( (line = br2.readLine()) != null) {
                         st.append(line+"\n");                            
                         System.out.println(line);
                     }
                  }
        } catch (Exception ex) {
            System.out.println(" error unable to execute the test... Error message follow:");
		ex.printStackTrace();            
        }
 	Writer out = new BufferedWriter(new OutputStreamWriter(System.out));
        //String[] args = {"/bin/bash", "-C", "executable/seqboot", ">executable/Inputs/inputSB", ">sortieSB"};

	
	
       
	
        File fwFile = new File(current+"/Result/Phylip");
        fwFile.mkdir();
        File fwFileT = new File(current+"/Result/Phylip/"+NameGene);
        fwFileT.mkdir();
        File fwFileCompt = new File(current+"/Result/Phylip/"+NameGene+"/"+pasWin);
        fwFileCompt.mkdir();
        Copy(false, NameGene, "outfile", "outfileSB", pasWin);
        File file = new File("outfile");
        file.delete();
        File fileSortie = new File("sortieSB");
        fileSortie.delete();
	
    }


    public void ProtDist(String NameGene, String pasWin) throws IOException, InterruptedException{
	try {
	StringBuilder st=new StringBuilder(); 
        Runtime runtime = Runtime.getRuntime();	
	File fi=new File("");
	createCommandFile("executeProtdist", "#!/bin/sh\nexecutable/protdist <"+fi.getAbsolutePath()+"/executable/Inputs/inputProtD >sortieProtD\n");
        runtime.exec("chmod +x executeProtdist");
	String[] command=new String[5];
        for (int i=0; i<command.length;i++) command[i]="";
        //Execute the Makefile_test.command.sh  
        File f=new File("");
        command[0]=f.getAbsolutePath()+"/executeProtdist";
              
	       Process p = runtime.exec(command);
               InputStream stdoutput = p.getInputStream();
               InputStream stderr = p.getErrorStream();
               InputStreamReader isr = new InputStreamReader(stdoutput);
               InputStreamReader isr2 = new InputStreamReader(stderr);
               BufferedReader br = new BufferedReader(isr);
               BufferedReader br2 = new BufferedReader(isr2);
               String line = null;
               //--System.out
               while ( (line = br.readLine()) != null) {
                 st.append(line+"\n");
                  System.out.println(line);
               }
                   int exitVal = p.waitFor();                                                                
                  //--System.err
                   if (exitVal!=0) {
                       System.out.println("Protdist error "+exitVal);
                      
                     while ( (line = br2.readLine()) != null) {
                         st.append(line+"\n");                            
                         System.out.println(line);
                     }
                  }
        } catch (Exception ex) {
            System.out.println(" error unable to execute the test... Error message follow:");
		ex.printStackTrace();            
        }
 	Writer out = new BufferedWriter(new OutputStreamWriter(System.out));



        Copy(false, NameGene, "outfile", "outfileProtD", pasWin);
        File file = new File("outfile");
        file.delete();
        File fileSortie = new File("sortieProtD");
        fileSortie.delete();

    }

    public void Neighbor(String NameGene, String pasWin) throws IOException, InterruptedException{

	try {
	StringBuilder st=new StringBuilder(); 
        Runtime runtime = Runtime.getRuntime();	
	File fi=new File("");
	createCommandFile("executeNeighbor", "#!/bin/sh\nexecutable/neighbor <"+fi.getAbsolutePath()+"/executable/Inputs/inputNJ >sortieNJ\n");
        runtime.exec("chmod +x executeNeighbor");
	String[] command=new String[5];
        for (int i=0; i<command.length;i++) command[i]="";
        //Execute the Makefile_test.command.sh  
        File f=new File("");
        command[0]=f.getAbsolutePath()+"/executeNeighbor";
              
	       Process p = runtime.exec(command);
               InputStream stdoutput = p.getInputStream();
               InputStream stderr = p.getErrorStream();
               InputStreamReader isr = new InputStreamReader(stdoutput);
               InputStreamReader isr2 = new InputStreamReader(stderr);
               BufferedReader br = new BufferedReader(isr);
               BufferedReader br2 = new BufferedReader(isr2);
               String line = null;
               //--System.out
               while ( (line = br.readLine()) != null) {
                 st.append(line+"\n");
                  System.out.println(line);
               }
                   int exitVal = p.waitFor();                                                                
                  //--System.err
                   if (exitVal!=0) {
                       System.out.println("Neighbor error "+exitVal);
                      
                     while ( (line = br2.readLine()) != null) {
                         st.append(line+"\n");                            
                         System.out.println(line);
                     }
                  }
        } catch (Exception ex) {
            System.out.println(" error unable to execute the test... Error message follow:");
		ex.printStackTrace();            
        }
 	Writer out = new BufferedWriter(new OutputStreamWriter(System.out));


        Copy(false, NameGene, "outtree", "outtreeNJ", pasWin);
        File file = new File("outfile");
        file.delete();
        File fileTree = new File("outtree");
        fileTree.delete();
        File fileSortie = new File("sortieNJ");
        fileSortie.delete();

    }


   public void Consence(String NameGene, String pasWin) throws IOException, InterruptedException{


	try {
	StringBuilder st=new StringBuilder(); 
        Runtime runtime = Runtime.getRuntime();	
	File fi=new File("");
	createCommandFile("executeConsense", "#!/bin/sh\nexecutable/consense <"+fi.getAbsolutePath()+"/executable/Inputs/inputCs >sortieCs\n");
        runtime.exec("chmod +x executeConsense");
	String[] command=new String[5];
        for (int i=0; i<command.length;i++) command[i]="";
        //Execute the Makefile_test.command.sh  
        File f=new File("");
        command[0]=f.getAbsolutePath()+"/executeConsense";
              
	       Process p = runtime.exec(command);
               InputStream stdoutput = p.getInputStream();
               InputStream stderr = p.getErrorStream();
               InputStreamReader isr = new InputStreamReader(stdoutput);
               InputStreamReader isr2 = new InputStreamReader(stderr);
               BufferedReader br = new BufferedReader(isr);
               BufferedReader br2 = new BufferedReader(isr2);
               String line = null;
               //--System.out
               while ( (line = br.readLine()) != null) {
                 st.append(line+"\n");
                  System.out.println(line);
               }
                   int exitVal = p.waitFor();                                                                
                  //--System.err
                   if (exitVal!=0) {
                       System.out.println("Consense error "+exitVal);
                      
                     while ( (line = br2.readLine()) != null) {
                         st.append(line+"\n");                            
                         System.out.println(line);
                     }
                  }
        } catch (Exception ex) {
            System.out.println(" error unable to execute the test... Error message follow:");
		ex.printStackTrace();            
        }
 	Writer out = new BufferedWriter(new OutputStreamWriter(System.out));


        Copy(false, NameGene, "outtree", "outtreeConsence", pasWin);
        File file = new File("outfile");
        file.delete();
        File fileSortie = new File("sortieCS");
        fileSortie.delete();
        File fileTree = new File("outtree");
        fileTree.delete();

    }

}
