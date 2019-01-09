/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.Date;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration; 

/**
 *
 * @author Nadia Tahiri
 */
public class Main {

    public static String [] listefichiers;
    public static String [] listeFiles;
    public static int nbFiles=0;
    public static int nbEspeces=0;
    public static int longueurSequence=0;
    public static int longueurMaxNomEspece=0;
    public static ArrayList<String> noms_Especes=new ArrayList<String>();
    public static String current=(new File("").getAbsolutePath());
    public  File repertoire=new File(current);
    public static File repertoireData=new File(current+"\\data_alignement");


    public static float moyBoostrapUtilisateur=0;
    public static float valSeuilRF=0;
    public static int tailleWindows=0;
    public static int pasWindows=0;
    public static String [] nomFileTree;
    public static int nb_SpeOutput =0;
    public static int nb_Spe [] = new int [2];

    /**
     * Pretraitement des differents fichiers d'entrees
     */
    public static void PretraitementEntrees() throws IOException {

         listerRepertoire(repertoireData);

         /**
          * Visualiser les noms des fichiers
          * de sorte qu'ils ne contiennent aucun espace
          * si certains le cas retourner la liste des noms des fichiers a corriger.
          */
         if (FileAvecEspace(listeFiles)){
             /**
              * si presence d'espace on va indiquer les noms des fichiers a modifier pour executer le progarmme sans bug
              */

             AfficheFileAvecEspace(listeFiles);
         }else{

             File fb = new File(current+"/Result");
             fb.mkdir();

             if (nbFiles<=1){
                 System.out.print("Votre dossier contient "+nbFiles+" fichier : \n");
             }else{
                 System.out.print("Votre dossier contient "+nbFiles+" fichiers : \n");
             }

             for (int j=0;j<nbFiles;j++){

                 File fwFile = new File(current+"/Result/"+listeFiles[j]);
                 fwFile.mkdir();

                 System.out.print(listeFiles[j]+"\n");

                  try{
                    System.out.print("---> Pretraitement terminé pour ce fichier.\n");
                }catch (Exception exception){
                    System.out.print ("---> Le fichier n'a pas été trouvé. "+exception+"\n");
                }
             }
         }

    }


    /**
     * vrai ou faux, le fichier existe ?
     */
     public static boolean fichierExiste(String nom_file) {

        boolean existe = false;
        try {
            File fichier = new File(nom_file);
            if (fichier.exists()) {
                existe = true;
            }
        } catch (Exception e) {
            existe = false;
        }
        return existe;

    }

     /*
      * Efface les noms des especes ainsi que leurs sequences
      * lorsqu'il y a des doublons des noms des espces
      * @params Nom du fichier (les sequences + nom especes en format fasta)
      * @output une autre fichier portant the same name mais sans l'extension (.txt) qu'en entrée.
      */
     public static void EffaceDoublon(String NameFile) throws IOException{

        FileWriter fw = new FileWriter("data_sequence/"+NameFile.substring(0, NameFile.length()-4));
        File f = new File ("data_sequence/"+NameFile);

        try{

            Scanner scanner = new Scanner (f);
            String nomEspeces;
            String sequence;

            int N=100;

            String []tmpNomEspeces=new String [N];
            for (int remp=0;remp<N;remp++){
                tmpNomEspeces[remp]="S.O.";
            }

            if (fichierExiste("data_sequence/"+NameFile)){
                int incrementation=0;
                while (scanner.hasNext()){
                    nomEspeces=scanner.next();
                    if (!AppartenanceTab(nomEspeces,tmpNomEspeces)){
                        tmpNomEspeces[incrementation]=nomEspeces;
                        fw.write(nomEspeces);
                        fw.write("\n");
                        sequence=scanner.next();
                        fw.write(sequence);
                        fw.write("\n");
                        incrementation++;
                    }else{
                        nomEspeces=scanner.next();
                    }

                    try{

                       }catch (NoSuchElementException exception){
                            break;
                    }
                }
                System.out.print("\nPretraitement terminé pour le fichier :\n");
            }
            scanner.close();
            fw.close();
        }catch (FileNotFoundException exception){
            System.out.println ("Le fichier n'a pas été trouvé");
        }

     }

    /**
     * Pour voir la presence d'une espece dans un tableau
     * ou plus generalement, pour voir la presence d'un element String dans un tableau
     */
     public static boolean AppartenanceTab(String nomEspeces,String [] tabNomEspeces) {

        boolean presence = false;
		
		int parc = 0;
        while (parc<tabNomEspeces.length && presence == false){
			if(tabNomEspeces[parc]!=null){
				if(tabNomEspeces[parc].equalsIgnoreCase(nomEspeces)){
					//Si on a rencontre l'espece dans le tableau
					//alors indiquer presence a TRUE
					presence=true;
				}
			}
			parc++;
        }
        return presence;

    }


    /**
     * Fonction qui va permettre de lister tous les fichiers dont l'extension est .txt sauf (Output.txt)
     * a partir du repertoire donne en entree
     */
    public static void listerRepertoire(File repertoire){

        int i;
        listefichiers=repertoire.list();
        listeFiles=new String[listefichiers.length];
        for(i=0;i<listefichiers.length;i++){
            if( !listefichiers[i].equalsIgnoreCase("Output.txt")||!listefichiers[i].equalsIgnoreCase("Tree")){
                listeFiles[nbFiles]=listefichiers[i];
                nbFiles++;
            }
        }

    }


    /**
     * Pour voir la presence d'un espace dans le nom des fichiers
     * @param
     */
    public static boolean FileAvecEspace(String [] listeFiles) {

        char space=' ';
        boolean espace = false;

        for (int lgtab=0;lgtab<nbFiles;lgtab++){
            espace = false;
            if (!listeFiles[lgtab].equalsIgnoreCase("null")){
                for (int lgstr=0;lgstr<listeFiles[lgtab].length();lgstr++){
                    if(listeFiles[lgtab].charAt(lgstr)==space){
                        espace=true;
                        return espace;
                    }
                }
            }else{
                lgtab=listeFiles.length+1;
            }
        }
        return espace;

    }


    /**
     * Indique les noms des fichiers qui contiennent un espace
     */
    public static void AfficheFileAvecEspace(String [] listeFiles) {

        char space=' ';
        for (int lgtab=0;lgtab<nbFiles;lgtab++){
            if (!listeFiles[lgtab].equalsIgnoreCase("null")){
                for (int lgstr=0;lgstr<listeFiles[lgtab].length();lgstr++){
                    if(listeFiles[lgtab].charAt(lgstr)==space){
                        System.out.print("Liste des noms de fichiers contenants des espaces :\n");
                        System.out.print(listeFiles[lgtab]+"\n");
                    }
                }
            }else{
                lgtab=listeFiles.length+1;
            }
        }

    }


    /**
     * Va permettre de compter le nombre d'especes et la longueur des alignements de sequences
     * Va permettre de determiner la taille du nom de l'espece le plus long
     * qui va permettre de structurer par la suite la structure du bon format
     */
    public static void ReformatterFile(File fasta) throws IOException{

         try{
             nbEspeces=0;
             Scanner sc = new Scanner (fasta);
             char sup='>';
             while (sc.hasNextLine()){

                 String ligne=sc.next();
                 char tm=ligne.charAt(0);
                 if (tm==sup){
                     if(longueurMaxNomEspece<ligne.length()){
                         longueurMaxNomEspece=ligne.length();
                     }
                    longueurSequence=0;
                    nbEspeces++;
                 }else{
                     longueurSequence+=ligne.length();
                 }
              }
              sc.close();
	  }catch (Exception exception){

          }

    }


    /**
     * Structure les fichiers d'entrees de PhylMl dans un bon format (phylip)
     */
    public static void InputPhyMl(File fasta, String fichierEntree) throws IOException{

        System.out.print(nbEspeces+"\t"+longueurSequence);
        FileWriter fastaFile = new FileWriter(new File(fichierEntree));
        fastaFile.write(nbEspeces+"\t"+longueurSequence);
        longueurMaxNomEspece=10;
         try{
             Scanner sc = new Scanner (fasta);
             char sup='>';
             while (sc.hasNextLine()){

                 String ligne=sc.next();
                 char tm=ligne.charAt(0);
                 if (tm==sup){
                     if (ligne.length()<longueurMaxNomEspece){
                        System.out.print("\n"+ligne.substring(1));
                        fastaFile.write("\n"+ligne.substring(1));

                        for (int deb=ligne.length();deb<=longueurMaxNomEspece;deb++){
                            System.out.print(" ");
                            fastaFile.write(" ");
                        }
                     }else{
                         if (ligne.length()==longueurMaxNomEspece){
                             System.out.print("\n"+ligne.substring(1)+" ");
                             fastaFile.write("\n"+ligne.substring(1)+" ");
                         }else{
                             System.out.print("\n"+ligne.substring(1, longueurMaxNomEspece+1));
                             fastaFile.write("\n"+ligne.substring(1, longueurMaxNomEspece+1));
                         }
                     }
                     System.out.print(" ");
                     fastaFile.write(" ");
                 }else{
                     System.out.print(ligne);
                     fastaFile.write(ligne);
                 }
              }
              sc.close();
              fastaFile.close();
          }catch (Exception exception){

          }
          fastaFile.close();

    }


    /**
     * Fonction qui va faire appel a PhyMl en ligne de commande
     * avec pour parametres:
     * en entrees: des séquences sequentielles sous format phylipp et il s'agit des acides amines
     * boostrap=100
     *///executable/PhyML_3.0_linux64 -i Result/rhodopsin/55_1/rhodopsin_55_1 -d aa -b 100 > tmp.txt
    public static void PhyMl(String fichier) throws IOException{

        Runtime runtime = Runtime.getRuntime();
        //Process process = runtime.exec("executable/PhyML_3.0_linux64 -i "+fichier+" -d aa -b 100 > tmp.txt");
		Process process = runtime.exec("executable/PhyML_3.0_linux64 -i "+fichier+" -d aa -b 100");
		System.out.println("Hello ");
        InputStream is = process.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;

         while ((line = br.readLine()) != null) {
			System.out.println(line);
        } 
    }

    /**
     * Va permettre de récuperer le nombre d'especes et la longueur des alignements de sequences
     * Va permettre de determiner la taille du nom de l'espece le plus long
     */
    public static void Comptage(File fichierEntree)throws IOException{

         try{
             Scanner sc = new Scanner (fichierEntree);
             nbEspeces=sc.nextInt();
             longueurSequence=sc.nextInt();
             int []longeurMaxEspece=new int[nbEspeces];
             int compteur=-1;
             longueurMaxNomEspece=0;

             //String []tabLigneResult;

             while (sc.hasNextLine()){
                 compteur++;
                 String ligne1=sc.next();
                 String ligne2=sc.next();

                  longeurMaxEspece[compteur]=ligne1.length();
                  if(longueurMaxNomEspece<longeurMaxEspece[compteur]){
                    longueurMaxNomEspece=longeurMaxEspece[compteur];
                  }
              }
              sc.close();
	  }catch (Exception exception){

          }

    }


     /**
     * copie le fichier source dans le dossier souhaite.
     * retourne vrai si cela réussit.
     */
    public static boolean Copy(boolean modeAjout,String source, String dest){

        try{
            FileWriter destFile = new FileWriter(dest,modeAjout);

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



    public static ArrayList<String> WindowsTrees(String nomFile,String fichierEntree, int tailleWindows,int pasWindows) throws IOException{


        ArrayList<String> noms_fichiers=new ArrayList<String>();
        try{
            File Input=new File (fichierEntree);
            Comptage(Input);

/*             System.out.print("longueurMaxNomEspece "+longueurMaxNomEspece+"\n");
            System.out.print("longueurSequence "+longueurSequence+"\n");
            System.out.print("nbEspeces "+nbEspeces+"\n"); */

            File fbR = new File(current+"/Result");
            fbR.mkdir();

            File fbF = new File(current+"/Result/"+nomFile);
            fbF.mkdir();

            File fbAll = new File(current+"/Result/"+nomFile+"/"+longueurSequence+"_1");
            fbAll.mkdir();

            //Copie le fichier de l'alignement total dans le repertoire adequate soit ("\\Result\\"+nomFile+"\\"+longueurSequence+"_1\\"+nomFile+"_"+longueurSequence+"_1").
            if (Copy(false, "data_alignement/"+nomFile, "Result/"+nomFile+"/"+longueurSequence+"_1/"+nomFile+"_"+longueurSequence+"_1")){
				/* System.out.print("La copy de /data_alignement/"+nomFile+" vers Result/"+nomFile+"/"+longueurSequence+"_1/"+nomFile+"_"+longueurSequence+"_1 a ette correctement effectue.\n"); */
                //Ajout du fichier avec l'alignement total dans notre liste de fichiers.
				noms_fichiers.add("Result/"+nomFile+"/"+longueurSequence+"_1/"+nomFile+"_"+longueurSequence+"_1");
            }/* else{
                System.out.print("Une erreur est survenue lors de la copy de /data_alignement/"+nomFile+" vers Result/"+nomFile+"/"+longueurSequence+"_1/"+nomFile+"_"+longueurSequence+"_1.\n");
            } */


            for (int i=longueurMaxNomEspece+1;i<=longueurSequence+longueurMaxNomEspece+1-tailleWindows;i+=pasWindows){

                    File fwFile = new File(current+"/Result/"+nomFile+"/"+tailleWindows+"_"+(i-longueurMaxNomEspece));
                    fwFile.mkdir();

                    noms_fichiers.add("Result/"+nomFile+"/"+tailleWindows+"_"+(i-longueurMaxNomEspece)+"/"+nomFile+"_"+tailleWindows+"_"+(i-longueurMaxNomEspece));
                    /* System.out.print("nomFile "+"Result/"+nomFile+"/"+tailleWindows+"_"+(i-longueurMaxNomEspece)+"/"+nomFile+"_"+tailleWindows+"_"+(i-longueurMaxNomEspece)+"\n"); */

                    FileWriter fastaFile = new FileWriter(current+"/Result/"+nomFile+"/"+tailleWindows+"_"+(i-longueurMaxNomEspece)+"/"+nomFile+"_"+tailleWindows+"_"+(i-longueurMaxNomEspece));
                    //Permet d'indiquer la bonne taille de l'alignement
                    if ((tailleWindows+i)<(longueurSequence+longueurMaxNomEspece+1)){
                        /* System.out.print(nbEspeces+"\t"+tailleWindows); */
                        fastaFile.write(nbEspeces+"\t"+tailleWindows);
                    }else{
                        /* System.out.print(nbEspeces+"\t"+(longueurSequence+longueurMaxNomEspece+1-i)); */
                        fastaFile.write(nbEspeces+"\t"+(longueurSequence+longueurMaxNomEspece+1-i));
                    }


                    //lecture du fichier texte
                    try{
                        InputStream ips=new FileInputStream(fichierEntree);
                        InputStreamReader ipsr=new InputStreamReader(ips);
                        BufferedReader br=new BufferedReader(ipsr);
                        String ligne;
                        ligne=br.readLine();

                        while ((ligne=br.readLine())!=null){
                            /* System.out.print("\n"); */
                            fastaFile.write("\n");
                            for (int deb=0;deb<longueurMaxNomEspece;deb++){
                                /* System.out.print(ligne.charAt(deb)); */
                                fastaFile.write(ligne.charAt(deb));
                            }

                            /* System.out.print(" "); */
                            fastaFile.write(" ");
                            if ((tailleWindows+i)<=(longueurSequence+longueurMaxNomEspece+1)){
                                for (int j=i;j<i+tailleWindows;j++){
                                    /* System.out.print(ligne.charAt(j)); */
                                    fastaFile.write(ligne.charAt(j));
                                }
                            }else{
                                int ind=i;
                                while (ind<(longueurSequence+longueurMaxNomEspece+1)){
                                    /* System.out.print(ligne.charAt(ind)); */
                                    fastaFile.write(ligne.charAt(ind));
                                    ind++;
                                }
                            }
                        }
                        /* System.out.print("\n\n"); */
                        br.close();
                    }

                    catch (Exception e){
                        System.out.println(e.toString());
                    }
                    fastaFile.close();

            }

            }catch (Exception exception){

            }
            return noms_fichiers;
    }



    /**
     * copie le fichier source dans le dossier Result
     * retourne vrai si cela réussit.
     */
    public static boolean copyFile(boolean modeAjout,String source, String dest){

        try{
            File fwFileT = new File(current+"/Result/Trees/"+dest.substring(0,dest.length()-4));
            fwFileT.mkdir();
            FileWriter destFile = new FileWriter(current+"/Result/Trees/"+dest.substring(0,dest.length()-4)+"/"+dest,modeAjout);

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

	/**
     * copie le fichier source dans le dossier Result
     * retourne vrai si cela réussit.
     */
    public static boolean copyFileTrees(boolean modeAjout,String source, String dest){

        try{
			File fwFileT = new File(current+"/Result/Trees/"+dest.substring(0,dest.length()-4));
            fwFileT.mkdir();
            FileWriter destFile = new FileWriter(current+"/Result/Trees/"+dest.substring(0,dest.length()-4)+"/"+dest,modeAjout);

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

    public static void AffichageResult (int n, String nomFile, int tailleWindows, int pasWindows, ArrayList<Integer> indicesConservees, float []moyBoostrapArbre, String output, int nbEspcesGeo,String paramClima, double []robinson_foulds, double [] bootstrap) throws FileNotFoundException, IOException{


	FileWriter TabResult = new FileWriter(current+"/Result/synthese_"+nomFile+".txt",true);
    String egal=     "============================================================";
    String petitTire="------------------------------------------------------------";

    Date maDate = new Date();
	// affichage:
	System.out.print("Date d'execution :"+maDate.toString()+"\n");
 	TabResult.write("Date d'execution :"+maDate.toString()+"\n");

	String enteteGene="+-------------------------------------------+\n"+
                      "|                    GENE                   |\n"+
                      "+-------------------------------------------+";

        System.out.print(enteteGene+"\n");
        TabResult.write(enteteGene+"\n");

        System.out.print("|Nom Gene : "+nomFile+"\n");
        TabResult.write("|Nom Gene : "+nomFile+"\n");

        System.out.print("|Longueur Gene : "+longueurSequence+"\n");
        TabResult.write("|Longueur Gene : "+longueurSequence+"\n");

        System.out.print("|Nombre Especes : "+nbEspeces+"\n");
        TabResult.write("|Nombre Especes : "+nbEspeces+"\n");

        System.out.print("\n");
        TabResult.write("\n");

        String enteteWindows="+-------------------------------------------+\n"+
                             "|                   WINDOWS                 |\n"+
                             "+-------------------------------------------+";
        System.out.print(enteteWindows+"\n");
        TabResult.write(enteteWindows+"\n");

        System.out.print("|Taille de la fenetre : "+tailleWindows+"\n");
        TabResult.write("|Taille de la fenetre : "+tailleWindows+"\n");

        System.out.print("|Pas d'avancement de la fenetre : "+pasWindows+"\n");
        TabResult.write("|Pas d'avancement de la fenetre : "+pasWindows+"\n");

        System.out.print("\n");
        TabResult.write("\n");

		String enteteTrees="+-------------------------------------------+\n"+
                           "|                   CLIMAT                   |\n"+
                           "+-------------------------------------------+";
        System.out.print(enteteTrees+"\n");
        TabResult.write(enteteTrees+"\n");

        System.out.print("|=> CLIMATIC TREE : "+"\n");
        TabResult.write("|=> CLIMATIC TREE : "+"\n");
		
		System.out.print("\t"+"Parametre climatique:"+paramClima+"\n");
        TabResult.write("\t"+"Parametre climatique:"+paramClima+"\n");

		System.out.print("\t"+"Nombre d'especes sur l'arbre de distribution:"+nbEspcesGeo+"\n");
        TabResult.write("\t"+"Nombre d'especes sur l'arbre de distribution:"+nbEspcesGeo+"\n");

        System.out.print("\t"+"Nombre d'especes identiques entre l'arbre de distribution et l'arbre d'especes:"+n +"\n");
        TabResult.write("\t"+"Nombre d'especes identiques entre l'arbre de distribution et l'arbre d'especes:"+n +"\n");
		

        System.out.print("\n");
        TabResult.write("\n");

		
        System.out.print("Result :\n");
        TabResult.write("Result :\n");

        System.out.print(egal+"\n");
        TabResult.write(egal+"\n");

        System.out.print("Positions"+"\t|");
        TabResult.write("Positions"+"\t|");

        System.out.print("RF Normalise"+"\t");
        TabResult.write("RF Normalise"+"\t");
		
        System.out.print("|");
        TabResult.write("|");

		System.out.print("Boostrap Moyen"+"\t");
        TabResult.write("Boostrap Moyen"+"\t");
		
        System.out.print("|");
        TabResult.write("|");

        System.out.print("\n");
        TabResult.write("\n");

        System.out.print(petitTire+"\n");
        TabResult.write(petitTire+"\n");
        String formatter;

        for (int deBis:indicesConservees){
            //Calcul des positions:
            int posDebut=1+((deBis-2)*pasWindows);
            int posFin=posDebut+tailleWindows-1;
            if (deBis==1){
                posDebut=1;
                posFin=longueurSequence;
                formatter = String.format("%1$5s %2$-1s %3$-5s","*"+posDebut, "-", ""+posFin);
                System.out.print(formatter+"\t");
                TabResult.write(formatter+"\t");
            }else{
                if(posFin<longueurSequence){
                   formatter = String.format("%1$5s %2$-1s %3$-5s",""+posDebut, "-", ""+posFin);
                   System.out.print(formatter+"\t");
                   TabResult.write(formatter+"\t");
                }else{
                   formatter = String.format("%1$5s %2$-1s %3$-5s",""+posDebut, "-", ""+longueurSequence);
                   System.out.print(formatter+"\t");
                   TabResult.write(formatter+"\t");
                }
            }
            System.out.print("|");
            TabResult.write("|");
		    
			double RFnormalise =robinson_foulds[deBis];
			formatter = String.format("%1$-21s",""+RFnormalise);
			System.out.print(formatter+"\t");
			TabResult.write(formatter+"\t");

            System.out.print("|");
            TabResult.write("|");

			formatter = String.format("%1$-23s",""+bootstrap[deBis]);
			System.out.print(formatter+"\t");
			TabResult.write(formatter+"\t");

            System.out.print("|");
            TabResult.write("|");
            System.out.print("\n");
            TabResult.write("\n");

        }
        System.out.print(egal+"\n\n");
        TabResult.write(egal+"\n\n");

        TabResult.close();
    }


     public static double RecuperationRF (String nomFile) throws FileNotFoundException, IOException{
        double RFnormalise = -1;
		double RF = 0;
		File fi = new File ("");
		int nbEspecesIdentiques = 0;
        String outFile=fi.getAbsolutePath()+"/Result/Trees/"+nomFile;
        //System.out.println(outFile+"\n");
        File out = new File (outFile);
        Scanner scanna= new Scanner(out);
        String ligResult=scanna.next();
        String []tabLigResult=new String[4];
        for (int k=0;k<1;k++){
            tabLigResult=ligResult.split("<>");
            RF = Double.parseDouble(tabLigResult[1]);	
			nbEspecesIdentiques = Integer.parseInt(tabLigResult[0]);
            if(nbEspecesIdentiques!=3){
                RFnormalise = (RF/(2*nbEspecesIdentiques-6));
				//RFnormalise = RF;
            }else{
                System.out.print("Erreur: Division par 0!\n");
            }
        }
        scanna.close();
		//System.out.println("La distance de Robinson-Foulds normalise est : "+RFnormalise+".\n");
        return RFnormalise;
    }

	public static double RecuperationBD (String nomFile) throws FileNotFoundException, IOException{
		double BD = 0;
		File fi = new File ("");
        String outFile=fi.getAbsolutePath()+"/Result/Trees/"+nomFile;
        //System.out.println(outFile+"\n");
        File out = new File (outFile);
        Scanner scanna= new Scanner(out);
        String ligResult=scanna.next();
        String []tabLigResult=new String[4];
        for (int k=0;k<1;k++){
            tabLigResult=ligResult.split("<>");
            BD = Double.parseDouble(tabLigResult[2]);	
        }
        scanna.close();
		//System.out.println("La distance de Robinson-Foulds normalise est : "+RFnormalise+".\n");
        return BD;
    }

     public static void Impression(int nb_trees, String [] nomFileTree, float moyBoostrapUtilisateur, float valSeuilRF, int tailleWindows, int pasWindows){


        NodeTree node = new NodeTree ();

        for (int fileTree=0; fileTree<nomFileTree.length;fileTree++){
        do{
            nb_SpeOutput=node.RecupererNomSpecies("arbres/"+nomFileTree[fileTree], (fileTree+1));
			nb_Spe[fileTree]=nb_SpeOutput;
            /* System.out.print("Votre arbre de distribution "+nomFileTree[fileTree]+" contient: "+nb_SpeOutput+" especes\n"); */
            /* System.out.print("Votre fichier Output_"+(fileTree+1)+" se trouvant au dossier tmp de votre projet a ete correctement cree et contient tous les noms des especes de votre arbre.\n"); */
         }while (!fichierExiste("arbres/"+nomFileTree[fileTree]));
            /* System.out.print("-----------------------------\n"); */
        }
    }




    public static int NbEpseceEnCommun ( int nb_SpeAlignement) throws FileNotFoundException{
		System.out.println("hello");
        int n = 0;  //nb d'especes identiques qux deux arbres
        //Remplissage des noms d'especes de l'arbre etudie.
        File f1 = new File ("tmp/Output_1");
        Scanner sca = new Scanner (f1);
		System.out.println("nb_SpeAlignement "+nb_SpeAlignement);
        String []nomEspecesTree = new String [nb_SpeAlignement];
        int cmp =0;

        while (sca.hasNext()){
            nomEspecesTree[cmp]=sca.next();
			System.out.println("cmp "+cmp);
            cmp ++;
        }
		System.out.println("hsgshs "+cmp);
        //Fermeture des fichiers
        sca.close();


        //Ouverture du fichier contenant la liste des noms d'especes de l'arbre geographiques.
        //String fichierDataListeEspeces = new String ("tmp\\Output_"+de);
        File f2 = new File ("tmp/Output_0");
        Scanner scannage = new Scanner (f2);

		String esp_en_cours = "";
        while (scannage.hasNext()){
			esp_en_cours = scannage.next();
            if(AppartenanceTab(esp_en_cours, nomEspecesTree)){
                n++;
            }
        }

        System.out.println("Le nombre d'especes identiques entre les deux arbres (Output_1 et Output_0): "+n+"\n");
        return n;
    }

    /*
     * Procedure du mode operatoire pour la réalisation des windows-Trees
     * suivit de la construction des arbres pour chaque fenetre
     * puis des calculs des differentes distances
     */
    public static void Run_Process (String nomFileGeo,String nomFileGene,int nbEspcesGeo,String paramClima) throws IOException, InterruptedException{

        File fi = new File("");
        Scanner kb = new Scanner(System.in);
        boolean auMoinsUnArbre = false;
        boolean auMoinsUnArbre_G = false;
        ArrayList<Integer> indicesConservees = new ArrayList<Integer>();//Conservation des indices que nous allons conserver.
        ArrayList<Integer> indices = new ArrayList<Integer>();//Tous les indices.
        ArrayList<String> noms_fichiers = new ArrayList<String>();
        HGT hgt = new HGT ();
        Boostrap boostrap = new Boostrap();
        NodeTree node = new NodeTree ();
        int nb_SpeAlignement = 0;
        int n = 0; //nb d'especes identiques aux deux arbres
        double rfConsence = 0;
		boolean trouve = true;
		boolean error = false;
		double rf = 0;
		

		//int n[] = new int[2];

        int nb_windows=0; //Cette variable va stocker des fenetres à traiter.

        //Recuperation des differents fichiers contenant les fenetres des alignements sous format phylipp dans un tableau nom_fichiers.
        noms_fichiers = WindowsTrees(nomFileGene,"data_alignement/"+nomFileGene,tailleWindows,pasWindows);
	
        int longueurSequenceMoinsPas=longueurSequence;//Variable qui va indiquer la taille de la nouvelle sous sequence
        while(tailleWindows<=longueurSequenceMoinsPas){
            nb_windows++;
            longueurSequenceMoinsPas-=pasWindows;
        }
		
		double robinson_foulds [] = new double [nb_windows+2];
		double bootstrap [] = new double [nb_windows+2];
		
        float []moyBoostrapArbre=new float[nb_windows+1];
        float []moyBoostrapArbreConsence=new float[nb_windows+1];

        //Creation du dossier Trees sous le repertoire Result.
        File fwFile = new File(current+"/Result/Trees");
        fwFile.mkdir();

        //Traitement pour l'ensemble des fenetres et plus un pour l'alignement au complet.
        for (int de=0;de<(nb_windows+1);de++){
            indices.add(de+1);
            auMoinsUnArbre_G = true;

            /**
             * PHYLIP
             */
            String pasWin=""+(de+1);
            //creation des input pour les executions des executables de Phylip.
            FilesInputPhylip fileInput = new FilesInputPhylip();
            fileInput.InputSB(noms_fichiers.get(de));
            fileInput.InputProtD(nomFileGene, pasWin);
            fileInput.InputNJ(nomFileGene, pasWin);
            fileInput.InputCs(nomFileGene, pasWin);

            //Execution des executable phylip
            Phylip phylip = new Phylip();
            phylip.SeqBoot(nomFileGene, pasWin);
            phylip.ProtDist(nomFileGene, pasWin);
            phylip.Neighbor(nomFileGene, pasWin);
            phylip.Consence(nomFileGene, pasWin);

            //Formatter correctement la sortie du Tree de Consence
            fileInput.FileCs(nomFileGene, pasWin);

            //appel de la fonction boostrap afin de calculer le boostrap moyen de l'arbre issu de Consence
            moyBoostrapArbreConsence[de]=boostrap.RecupererMoyenneBoostrapConsence("Result/Phylip/"+nomFileGene+"/"+pasWin+"/outtreeConsence");

			System.out.println(nomFileTree.length);
            //for (int fileTree=0; fileTree<nomFileTree.length; fileTree++){
                //Copie du premier arbre Tree du parametre d'entree, afin de preparer le fichier d'entree de la fonction Criterion
				nb_SpeAlignement = node.RecupererNomSpecies("arbres/"+nomFileGeo, 0);
				copyFile(false,"arbres/"+nomFileGeo, nomFileGene+"_"+nomFileGeo+".txt");

                //Copie du deuxieme arbre issu de l'execution de PhyMl, afin de preparer le fichier d'entree de la fonction Criterion
                nb_SpeAlignement = node.RecupererNomSpecies("Result/Phylip/"+nomFileGene+"/"+pasWin+"/outtreeConsences", 1);
				copyFile(true,"Result/Phylip/"+nomFileGene+"/"+pasWin+"/outtreeConsences", nomFileGene+"_"+nomFileGeo+".txt");

                //Bon emplacement de réaliser l'execution de RF ou HGT
                //System.out.print("File PPPPP "+nomFile+"Cs_"+(de+1)+"_"+(fileTree+1)+".txt\n");
                error = hgt.CriterionConsence(nomFileGene+"_"+nomFileGeo+".txt");
				
                //Revoir l'erreur
                //MESSAGE:SEGMENTATION FAULT #11 DETECTED
                //Use valgrind or gdb to fix the problem

                //Appel de la fonction qui a pour entree deux fichiers, contenant les noms des especes (tmp\\Output_ )
                //et retourne le nombre d'especes en commun des 2 arbres (variable n)
                //n [fileTree]= NbEpseceEnCommun(nb_SpeAlignement, fileTree); //nb d'especes identiques aux deux arbres
				if(de==0){ //pour calculer qu'une fois n
					n = NbEpseceEnCommun(nb_SpeAlignement); //nb d'especes identiques aux deux arbres
					System.out.println("N ===  "+n);
				}
                if (error){
                    rfConsence = valSeuilRF;
                    String outFile=fi.getAbsolutePath()+"/Result/Trees/"+nomFileGene+"_"+nomFileGeo+"/outputConsence.txt";
                    FileWriter file = new FileWriter(outFile);
                    file.write(""+rfConsence);
                    file.close();
                }else{
                    rfConsence = RecuperationRF (nomFileGene+"_"+nomFileGeo+"/outputConsence.txt");
					robinson_foulds[de+1]=rfConsence;
                }

            //}

			//PhyMl(noms_fichiers.get(de));

			 moyBoostrapArbre[de]=boostrap.RecupererMoyenneBoostrap(noms_fichiers.get(de)+"_phyml_tree.txt");
			 bootstrap[de+1]=moyBoostrapArbre[de];
			 //if (moyBoostrapArbre[de]>=moyBoostrapUtilisateur){
			 if (moyBoostrapArbre[de]>=moyBoostrapUtilisateur){
				//pour eviter de mettre plusieurs foisles mêmes postions dans le tableau indicesConservees.
/* 				if (trouve){
					indicesConservees.add(de+1);
					trouve = false;
				} */
				indicesConservees.add(de+1);
				auMoinsUnArbre=true;
				//Bon emplacement pour réaliser l'execution de RF ou de HGT
				hgt.Criterion(nomFileGene+"_"+nomFileGeo+".txt");
			
				//Appel de la fonction qui a pour entree deux fichiers, contenant les noms des especes (tmp\\Output_ )
				//et retourne le nombre d'especes en commun des 2 arbres (variable n)

				rf = RecuperationRF (nomFileGene+"_"+nomFileGeo+"/output.txt");

			}
        }

        if(auMoinsUnArbre){
            //Synthese des resultats valides
            AffichageResult( n, nomFileGene, tailleWindows, pasWindows, indicesConservees, moyBoostrapArbre, "output.txt",nbEspcesGeo,paramClima,robinson_foulds,bootstrap);
        }


        if(auMoinsUnArbre_G){
           //Synthese de tous les resultats
            //AffichageResult( n, nomFile+"Cs", tailleWindows, pasWindows, indices, moyBoostrapArbreConsence, "outputConsence.txt",nbEspcesGeo,paramClima);
        }

    }


	public static double Run_Distance(String nomFileGen, String nomFileGeo, int choix) throws IOException, InterruptedException{
			HGT hgt = new HGT ();
			double distance = 0;
			//Copie du premier arbre Tree du parametre d'entree, afin de preparer le fichier d'entree de la fonction Criterion
			if (!copyFile(false,"arbres/"+nomFileGen, nomFileGen+"_"+nomFileGeo+".txt")){
				distance = -2;
			}

			//Copie du deuxieme arbre issu de l'execution de PhyMl, afin de preparer le fichier d'entree de la fonction Criterion
			if (!copyFile(true,"arbres/"+nomFileGeo, nomFileGen+"_"+nomFileGeo+".txt")){
				distance = -2;
			}
			//Bon emplacement pour réaliser l'execution de RF ou de HGT
			hgt.Criterion(nomFileGen+"_"+nomFileGeo+".txt");

			if(choix==1){
				distance = RecuperationRF (nomFileGen+"_"+nomFileGeo+"/output.txt");
			}else if(choix==2){
				distance = RecuperationBD (nomFileGen+"_"+nomFileGeo+"/output.txt");
			}
			
			return distance;
	}
	
	public static double [][] afficher_matrice(int nb_trees, String [] nomFileTree, int nb_alignemnts, String [] nomFile,  int choix) throws IOException, InterruptedException{
		//stocker les resultats de distances dans une matrice a 2 dimensions (nb_alignemnts x nb_trees)
		double [][] mat_distance = new double [nb_alignemnts][nb_trees];
		
		System.out.print(" \t");
		for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
			System.out.print(nomFileTree[fileGeo]+"\t");
		}
		System.out.print("\n");
		
		double distance = 0;
		
		for (int fileGen=0;fileGen<nb_alignemnts;fileGen++){
			System.out.print(nomFile[fileGen]+"\t");
			for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
				distance = Run_Distance(nomFile[fileGen],nomFileTree[fileGeo],choix);
				mat_distance[fileGen][fileGeo] = distance;
				System.out.print(distance+"\t");
			}
			System.out.print("\n");
        }
		
		return mat_distance;
	}
	
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException {
        // TODO code application logic here

        //Lecture des parametres d'entrees, depuis un fichier.
        File parametres = new File(args[0]);
        Scanner scan = new Scanner(parametres);


        //Data of the File in input
        String []nomFile;
        int nb_alignemnts=scan.nextInt();

        nomFile=new String[nb_alignemnts];
        for (int file=0; file<nomFile.length;file++){
            do{
                nomFile [file] = scan.next();
				
				if (!fichierExiste("arbres/"+nomFile[file])){
                    System.out.println("arbres/"+nomFile[file]+" : Votre fichier genetique n'ai pas trouvable!\n");
					nomFile [file] = null;
                }

            }while(!fichierExiste("arbres/"+nomFile[file]));
        }

        int nb_trees = scan.nextInt();
        nomFileTree=new String[nb_trees];
        for (int fileT=0; fileT<nomFileTree.length;fileT++){
            do{
                nomFileTree [fileT] = scan.next();
				if (!fichierExiste("arbres/"+nomFileTree[fileT])){
                    System.out.println("arbres/"+nomFileTree[fileT]+" : Votre fichier geographique n'ai pas trouvable!\n");
					nomFileTree [fileT] = null;
                }

            }while(!fichierExiste("arbres/"+nomFileTree[fileT]));
        }
		
		//parameters
		moyBoostrapUtilisateur = scan.nextFloat();
        valSeuilRF = scan.nextFloat();
        tailleWindows = scan.nextInt();
        pasWindows = scan.nextInt();
		
		System.out.println("\n\nSTATISTIQUES\n");
		//Pour calculer deux a deux les arbres génétiques et les arbres géographiques
		//Et affichage de la matrice RF
		
		//stocker les resultats des distances (RF normalisées, bipartition et quartet) dans une matrice a 2 dimensions (nb_alignemnts x nb_trees)
		double [][] mat_distance = new double [nb_alignemnts][nb_trees];
		
		Vector<Vector<String>> v_climatique = new Vector<Vector<String>>();
		Vector<String> v_climatique_unique = new Vector<String>();
		
		//vectors des minimums pour pour les proteines
		Vector<Vector<String>> v_gene = new Vector<Vector<String>>();
		Vector<String> v_gene_unique = new Vector<String>();
		
		System.out.println("\n\nMATRICE DES DISTANCES RF NORMALISEES");
		mat_distance = afficher_matrice(nb_trees,nomFileTree,nb_alignemnts,nomFile,1); 


		System.out.println("\n\nVECTEUR CONDITION CLIMATIQUE => LISTE DES GENES SIMILAIRES");
		//Connaitre les minimums pour chaque condition climatiques
		double minimum = 100000.0;
		
		for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
			minimum = 100000.0;
			//connaitre le minimum pour chaque condition climatique
			for (int fileGen=0;fileGen<nb_alignemnts;fileGen++){
				if(mat_distance[fileGen][fileGeo]<minimum){
					minimum = mat_distance[fileGen][fileGeo];
				}
			}
			
			//stocker tous les genes ayant ce minimum dans un vecteur
			v_climatique.add(new Vector<String>());
			for (int fileGen=0;fileGen<nb_alignemnts;fileGen++){				
				if(mat_distance[fileGen][fileGeo]==minimum){
					v_climatique.get(fileGeo).add(nomFile[fileGen]);
				}
			}			
		}
		
		//imprimer le vecteur
		for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
			System.out.println ((fileGeo+1)+") "+nomFileTree[fileGeo]+" => "+v_climatique.get(fileGeo));
		}
		
		int nbEspcesGeo = 0;
		//Appel de PhyML pour l'arbre geo nomFileTree[fileGeo] et l'arbre genetique v_climatique.get(fileGeo).get(fileGenetic).
		for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
			nbEspcesGeo = NodeTree.nbSpeciesArbres("arbres/"+nomFileTree[fileGeo]);
			for (int fileGenetic=0;fileGenetic<v_climatique.get(fileGeo).size();fileGenetic++){
				Run_Process(nomFileTree[fileGeo],v_climatique.get(fileGeo).get(fileGenetic),nbEspcesGeo,nomFileTree[fileGeo]);
			}
		}

		
		NodeTree node = new NodeTree ();
		//stocker les genes uniques dans un vecteur pour appeller PhyML
		for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
			nbEspcesGeo = NodeTree.nbSpeciesArbres("arbres/"+nomFileTree[fileGeo]);
			for (int fileGen=0;fileGen<v_climatique.get(fileGeo).size();fileGen++){
				if(!v_climatique_unique.contains(v_climatique.get(fileGeo).get(fileGen))){
					v_climatique_unique.add(v_climatique.get(fileGeo).get(fileGen));
					/* Run_Process(nomFileTree[fileGeo],v_climatique.get(fileGeo).get(fileGen),nbEspcesGeo,nomFileTree[fileGeo]); */
				}
			}
		}
		 
		System.out.println("\nLa liste des genes interessants:");
		System.out.println(v_climatique_unique);
		 
 		//Impression(nb_trees, nomFileTree, moyBoostrapUtilisateur, valSeuilRF, tailleWindows, pasWindows);
		//appel de PhyML
		System.out.println(v_climatique_unique.size());
 		for (int el=0;el<v_climatique_unique.size();el++){
			System.out.println(v_climatique_unique.get(el));
			//Run_Process (v_climatique_unique.get(el));
		}  
		
		// Remove all
		v_climatique.clear();
		
		
		
		
		System.out.println("\n\nVECTEUR DES PROTEINES => LISTE DES CONDITIONS CLIMATIQUES");
		//Connaitre les minimums pour chaque gene
		
		for (int fileGen=0;fileGen<nb_alignemnts;fileGen++){
			minimum = 100000.0;
			//connaitre le minimum pour chaque condition climatique
			for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
				if(mat_distance[fileGen][fileGeo]<minimum){
					minimum = mat_distance[fileGen][fileGeo];
				}
			}
			
			//stocker tous les genes ayant ce minimum dans un vecteur
			v_gene.add(new Vector<String>());
			for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){				
				if(mat_distance[fileGen][fileGeo]==minimum){
					v_gene.get(fileGen).add(nomFileTree[fileGeo]);
				}
			}			
		}
		
		//imprimer la vecteur
		for (int fileGen=0;fileGen<nb_alignemnts;fileGen++){
			System.out.println ((fileGen+1)+") "+nomFile[fileGen]+" => "+v_gene.get(fileGen));
		}
		
		//stocker les conditions climatiques uniques dans un vecteur pour appeller PhyML
		for (int fileGen=0;fileGen<nb_alignemnts;fileGen++){
			for (int fileGeo=0;fileGeo<v_gene.get(fileGen).size();fileGeo++){
				if(!v_gene_unique.contains(v_gene.get(fileGen).get(fileGeo))){
					v_gene_unique.add(v_gene.get(fileGen).get(fileGeo));
					
				}
			}
		}
		
		System.out.println(" ");
		System.out.println("La liste des conditions climatiques interessants:");
		System.out.println(v_gene_unique);
	
		// Remove all
		v_gene.clear();
				
		
		//PhyML pour les proteines interessantes
		/* Impression(nb_trees, nomFileTree, moyBoostrapUtilisateur, valSeuilRF, tailleWindows, pasWindows);
         for (int file=0;file<nb_alignemnts;file++){
            Run_Process (nomFile[file]);
        } */
		
		// Affichage simple du contenu à l'aide d'une boucle
/* 		for(int i=0; i<v_climatique.size(); i++){
			System.out.println("taille 2 = "+v_climatique.get(i).size());
			for(int j=0; j<v_climatique.get(i).size(); j++){
				System.out.println (i+" => "+v_climatique.get(i).get(j));
			}
		} */
 
		
		
		//Pour calculer deux a deux les arbres génétiques et les arbres géographiques
		//Et affichage de la matrice BD
		
		System.out.println("\n\nMATRICE DES DISTANCES DE BIDARTITIONS");
		//mat_distance = afficher_matrice(nb_trees,nomFileTree,nb_alignemnts,nomFile,2); 
		
		/* for (int fileGen=0;fileGen<nb_alignemnts;fileGen++){
			for (int fileGeo=0;fileGeo<nb_trees;fileGeo++){
				distance = Run_Distance(nomFile[fileGen],nomFileTree[fileGeo],choix);
				mat_RF_norm [fileGen][fileGeo] = distance;
			}
        } */
		 
        /* moyBoostrapUtilisateur = scan.nextFloat();
        valSeuilRF = scan.nextFloat();
        tailleWindows = scan.nextInt();
        pasWindows = scan.nextInt();

		
        Impression(nb_trees, nomFileTree, moyBoostrapUtilisateur, valSeuilRF, tailleWindows, pasWindows);
        for (int file=0;file<nb_alignemnts;file++){
            Run_Process (nomFile[file]);
        } */
        scan.close();
		
		System.out.println("FIN NORMALE DU PROGRAMME");
    }


}
