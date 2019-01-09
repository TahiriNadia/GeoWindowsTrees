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

/**
 *
 * @author nadia
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

       public void SeqBoot(String NameGene, String pasWin) throws IOException, InterruptedException{

        Runtime runtime = Runtime.getRuntime();
        Writer out = new BufferedWriter(new OutputStreamWriter(System.out));
        String[] args = { "executable/seqboot", "<executable/inputSB", ">sortieSB"};
        final Process process = runtime.exec(args);
        InputStream is = process.getInputStream();

        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;

        while ((line = br.readLine()) != null) {
         System.out.println(line);
        }

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
        
        Runtime runtime = Runtime.getRuntime();
        Writer out = new BufferedWriter(new OutputStreamWriter(System.out));
        String[] args = { "executable/protdist", "<executable/inputProtD", ">sortieProtD"};
        final Process process = runtime.exec(args);
        InputStream is = process.getInputStream();

        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;

        while ((line = br.readLine()) != null) {
         System.out.println(line);
        }

        Copy(false, NameGene, "outfile", "outfileProtD", pasWin);
        File file = new File("outfile");
        file.delete();
        File fileSortie = new File("sortieProtD");
        fileSortie.delete();
    }

    public void Neighbor(String NameGene, String pasWin) throws IOException, InterruptedException{

        Runtime runtime = Runtime.getRuntime();
        Writer out = new BufferedWriter(new OutputStreamWriter(System.out));
        String[] args = { "executable/neighbor", "<executable/inputNJ", ">sortieNJ"};
        final Process process = runtime.exec(args);
        InputStream is = process.getInputStream();

        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;

        while ((line = br.readLine()) != null) {
         System.out.println(line);
        }

        Copy(false, NameGene, "outtree", "outtreeNJ", pasWin);
        File file = new File("outfile");
        file.delete();
        File fileTree = new File("outtree");
        fileTree.delete();
        File fileSortie = new File("sortieNJ");
        fileSortie.delete();
    }


   public void Consence(String NameGene, String pasWin) throws IOException, InterruptedException{

        Runtime runtime = Runtime.getRuntime();
        Writer out = new BufferedWriter(new OutputStreamWriter(System.out));
        String[] args = { "executable/consense", "<executable/inputCs", ">sortieCS"};
        final Process process = runtime.exec(args);
        InputStream is = process.getInputStream();

        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;

        while ((line = br.readLine()) != null) {
         System.out.println(line);
        }

        Copy(false, NameGene, "outtree", "outtreeConsence", pasWin);
        File file = new File("outfile");
        file.delete();
        File fileSortie = new File("sortieCS");
        fileSortie.delete();
        File fileTree = new File("outtree");
        fileTree.delete();
    }

}
