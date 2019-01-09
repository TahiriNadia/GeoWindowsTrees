/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author nadia
 */
public class NodeTree {

    /*
     * Fonction permettant de recupere tous les noms des especes
     * a partir d'un arbre sous format newick
     * et retourne le nombre des especes presents
     */
    public int RecupererNomSpecies(String nomFileTree, int fileTree){
        int nb_spe=0;
        try{
                /**
                 * Pour vider le contenu du ficfier Output.txt au prealable.
                 */
                FileWriter outFile = new FileWriter("tmp/Output_"+fileTree);
                outFile.write("");
                outFile.close();
                File fichier= new File(nomFileTree);
                outFile = new FileWriter("tmp/Output_"+fileTree,true);
                 try {
                   FileReader flotLecture = new FileReader(fichier);
                   long longueurFichier= fichier.length();
                   int dejaLu = 0;
                   char car=0;
                   String nomSpecies="";


                   while (dejaLu < longueurFichier) {
                     car= (char)flotLecture.read();
                     dejaLu = dejaLu + 1;
                     if(!(car=='(' || car==')' || car==',' || car==':' || car==';'|| car=='.'|| car=='0' || car=='1' || car=='2' || car=='3' || car=='4' || car=='5' || car=='6' || car=='7' || car=='8' || car=='9')){
                         if (car=='_'){
                             nb_spe++;
                         }
                         nomSpecies+=car;
                     }else{
                         if((car==':')){
                             if(!nomSpecies.equalsIgnoreCase("")){
                                outFile.write(nomSpecies);
                                outFile.write("\n");
                                nomSpecies="";
                             }
                         }
                     }
                   }
                   outFile.close();
                   flotLecture.close();
                   return nb_spe;
                 } catch (IOException e) {
                 }

            }catch (Exception e){
            }
        return nb_spe;

    }

}
