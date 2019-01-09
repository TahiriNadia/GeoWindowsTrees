/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

/**
 *
 * @author nadia
 */
public class FilesInputPhylip {

     FileWriter fileInputSB = new FileWriter("executable/inputSB");
     FileWriter fileInputProtD = new FileWriter("executable/inputProtD");
     FileWriter fileInputNJ = new FileWriter("executable/inputNJ");
     FileWriter fileInputCs = new FileWriter("executable/inputCs");
    
    public FilesInputPhylip() throws IOException{

    }
    
    public void InputSB (String nomFolder) throws IOException{

        fileInputSB.write(nomFolder);
        fileInputSB.write("\nY");
        fileInputSB.write("\n73");
        fileInputSB.close();

    }
      
    public void InputProtD (String nomFolder, String pasWin) throws IOException{
        
        fileInputProtD.write("Result/Phylip/"+nomFolder+"/"+pasWin+"/outfileSB");
        fileInputProtD.write("\nP");
        fileInputProtD.write("\nP");
        fileInputProtD.write("\nP");
        fileInputProtD.write("\nM");
        fileInputProtD.write("\nD");
        fileInputProtD.write("\n100");
        fileInputProtD.write("\nY");
        fileInputProtD.close();
        
    }
    
    public void InputNJ (String nomFolder, String pasWin) throws IOException{
        
        fileInputNJ.write("Result/Phylip/"+nomFolder+"/"+pasWin+"/outfileProtD");
        fileInputNJ.write("\nM");
        fileInputNJ.write("\n100");
        fileInputNJ.write("\n73");
        fileInputNJ.write("\nY");
        fileInputNJ.close();
        
    }

    public void InputCs (String nomFolder, String pasWin) throws IOException{

        fileInputCs.write("Result/Phylip/"+nomFolder+"/"+pasWin+"/outtreeNJ");
        fileInputCs.write("\nY");
        fileInputCs.close();

    }

    public void FileCs (String nomFolder, String pasWin) throws IOException{

        File fileConsence = new File ("Result/Phylip/"+nomFolder+"/"+pasWin+"/outtreeConsence");
        FileWriter fileConsen = new FileWriter("Result/Phylip/"+nomFolder+"/"+pasWin+"/outtreeConsences");

        Scanner scan = new Scanner(fileConsence);

        int j=0;
        while (scan.hasNext()){
            String ligne=scan.next();
            fileConsen.write(ligne);
            //System.out.println(j+" "+ligne);
            j++;
        }

        fileConsen.close();
        fileConsence.exists();
    }
}
