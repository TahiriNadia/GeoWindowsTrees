/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author nadia
 */
public class Boostrap {

    /**
     * Fonction permettant de recupere tous les noms des especes
     * a partir d'un arbre sous format newick
     */
    public float RecupererMoyenneBoostrap(String nomFileTree){
        float moyBoostrap=0;
        try {
            File fichier= new File(nomFileTree);
            FileReader flotLecture = new FileReader(fichier);
            long longueurFichier= fichier.length();
            int dejaLu = 0;
            char car=0;

            String boostrap = "";

            int car_b=0;//Valeur binaire qui devient egal 1 en rencontrant une parenthese fermante
            int i_b=0;//compte le nombre de valeur de boostrap

            int sommeBoostrap=0;

            while (dejaLu < longueurFichier) {
                car= (char)flotLecture.read();
                dejaLu = dejaLu + 1;

                if (car==')'){
                    car_b=1;
                }else{
                    if(car!=':' && car_b==1){
                        boostrap= boostrap + car;
                    }
                }

                if(car==':' && car_b==1){
                    car_b=0;
                    sommeBoostrap+=Integer.parseInt(boostrap);
                    boostrap="";
                    i_b++;
                }

            }
            moyBoostrap=sommeBoostrap/i_b;
            System.out.println("Moyenne boostrap = "+moyBoostrap+"\n");
            flotLecture.close();
            return moyBoostrap;
        } catch (IOException e) {
        }
        return moyBoostrap;

    }

    /**
     * Fonction permettant de recupere tous les noms des especes
     * a partir d'un arbre sous format newick
     */
    public float RecupererMoyenneBoostrapConsence(String nomFileTree){
          float moyBoostrap=0;
        try {
            File fichier= new File(nomFileTree);
            FileReader flotLecture = new FileReader(fichier);
            long longueurFichier= fichier.length();
            int dejaLu = 0;
            char car;

            String boostrap = "";

            int car_b=0;//Valeur binaire qui devient egal 1 en rencontrant une parenthese fermante )
            int car_bp=0;//Valeur binaire qui devient egal 1 en rencontrant une parenthese fermante suivit de deux points ):
            int i_b=0;//compte le nombre de valeur de boostrap

            int sommeBoostrap=0;

            while (dejaLu < longueurFichier) {
                car= (char)flotLecture.read();
                dejaLu = dejaLu + 1;

                if (car==')'){
                    car_b=1;
                }

                if(car==':' && car_b==1){
                    car_bp=1;
                }

                if(car!=':' && car!='.' && car_bp==1){
                     boostrap= boostrap + car;
                }


                if(car=='.' && car_bp==1){
                    car_b=0;
                    car_bp=0;
                    sommeBoostrap+=Integer.parseInt(boostrap);
                    boostrap="";
                    i_b++;
                }

            }
            moyBoostrap=(sommeBoostrap/i_b);
            System.out.println("Moyenne boostrap = "+moyBoostrap+"\n");
            flotLecture.close();
            return moyBoostrap;
        } catch (IOException e) {
        }
        return moyBoostrap;

    }

}
