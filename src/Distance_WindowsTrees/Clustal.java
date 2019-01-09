/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 *
 * @author nadia
 */
public class Clustal {

     /**
     * Fonction qui va faire appel a ClustalW2 en ligne de commande
     * avec pour parametre:
     * une sortie en format FASTA
     */
    public void ClustalW(String fichier) throws IOException{

       Runtime runtime = Runtime.getRuntime();
       Process process = runtime.exec(new String[] {"executable\\clustalw2.exe", fichier,"-OUTPUT=FASTA"});
       InputStream is = process.getInputStream();
       InputStreamReader isr = new InputStreamReader(is);
       BufferedReader br = new BufferedReader(isr);
       String line;

       while ((line = br.readLine()) != null) {
         System.out.println(line);
       }
    }


}
