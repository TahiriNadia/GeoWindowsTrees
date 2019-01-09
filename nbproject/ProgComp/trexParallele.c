/*Approximation of the given dissimilarity DI by a tree metric DA*/
/*Copyright (C) 2004 by V. Makarenkov*/
/*Please see <http://alize.ere.umontreal.ca/~casgrain/en/labo/t-rex/> for more details*/

/*Enter the dissimilarity matrix's file name. Ther input file must be a text (ASCII) file,
containing the dissimilarity matrix in the PHYLIP format.
Values in the matrix should be separated by blanks (spaces) or tabs, as seen
in the sample files. Select one of the five reconstruction methods proposed. 

If you select the MW method you should also specify whether to use local
or global optimization, and define the matrix of weights. This latter matrix
can either be a file (same format as above) or a function of the dissimilarity
matrix: W = 1/(D^p) where W is the weight, D is the dissimilarity and
you specify p. If you enter 0.0 for p, then all weights W will be 1. 
*/

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "phylip.h"
#include "global.h"
#include "CalculDistances.c"
#include "DistanceMethode.c"
#include "seqboot.c"

/*Import de la librairie de MPI */
#include "mpi.h"

/*Definition de certaines constantes */
#define N1 1000
#define DIRDELIM '/' 
#define EOS '\0'

/*Variables de BioNJ*/
#define PREC 8  /* precision of branch-lengths  */
#define PRC  100
#define LEN  1000  /* length of taxon names */


/*Definition des structures*/
typedef struct word
{
	char name[LEN];
	struct word *suiv;
} WORD1;

typedef struct pointers
{
	WORD1 *head;
	WORD1 *tail;
} POINTERS;
/* End of BioNJ variables*/

/*variables passées en paramètre*/
int p_optionTR = 1;
int p_optionMiss = 2;
int p_method = 2;
int p_optionFunction = 1;
int p_iternumber = 5;
int p_option = 1;
int p_option1 = 1;
int p_p = 0;
int p_k = 5;
int p_boot = 2; /* 2 = no bootstrap , 0 = bootstrap , 1 = jackknife*/
int p_nbRep = 100;
char p_modele[100];
char p_input[100];
char p_output[100];
char p_stat[100];
char p_matrix[100];
char p_outboot[100];
char p_newick[100];
int m_sequenceMethod = 32887;	
double m_gapPenalty=0;  
int m_PEMV = 0;	
double m_paramA = -1;
int m_dataType = 1;   /*1 = distance , 2 = sequences*/
int totota;


/*Quelques variables globales*/ 
int n;
FILE *Output1;
FILE *TreeFile;
float method;
char * repertoire;
int kt;

/*Variables reliées à MPI*/
int    myid, numprocs, maxprocs;
int    namelen;
char   processor_name[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
double totaltime=0.0;


void construction(double**,double**,int*,double**);
void odp(double**,int*,int*,int*);
void parcour1(double**,double**,double**,int*,long int *, double *, int);
void parcour2(double**,double**,double**,int*,long int *, double *, int);
void NJ(double **,double **, int);
void ADDTREE(double **,double **, int);
void UNJ(double **,double **, int);
void construction2(double**,double**,int*);
int comparb(double**,double**);
int comparb1(double**,double**);
void coder(double**,int *,double*);
void decoder(double**,int *,double*);
void compute_criteres(double **,double **,double *,double **, int, char **);
void Tree_edges(double **,long int *,double *);
int floor1(double);
void approx_arb2(double**,double**,double**,double**, int*, long int*, double*);
void Scaling(double **, double *);
void Scaling1(double **,double **,double *,double *,double *);
void PrintEdges(long int *,double *, int, int, double*, double*, int);

float Agglomerative_criterion(int i, int j, float **delta, int r);
void Best_pair(float **delta, int r, int *a, int *b, int n);
float Branch_length(int a, int b, float **delta, int r);
void Compute_sums_Sx(float **delta, int n);
int Emptied(int i, float **delta);
void Finish(float **delta, int n, POINTERS *trees, FILE *output, double *l1, double *l2, double *l3, int *last1, int *last2, int *last3);
float Lamda(int a, int b, float vab, float **delta, int n, int r);
float Reduction10(int a, int b, int i, float lamda, float vab,float **delta);
float Reduction4(int a, float la, int b, float lb, int i, float lamda,float **delta);
float Variance(int i, int j, float **delta);
void BioNJ_main (double **DI, double **DA);

void RETICULATIONS (int n, double **DISS, double **D, long int *ARETE, double *LONGUEUR, int OptionFunction, double **W, int Iternumber, int *ReticulationsNumber, double *CRITERION1, double *CRITERION2, int);
int f_CalculerDistances(FILE *pfile,char *outFile, char **pmatNomsEspeces, double **pmatDistances,char **pmatSites,int pintNbreEspeces,int pintNbreSites,int pintIdMethod,double pdblPenalty,int pintPEMV,double pdblParamA,int ecrireResultats,int pblnWarningDejaAffiche);
int tryout(double **DI, double **W);
void LireOptions(char **argv, int argc);
void get_path(char *str,char *path);

/* Variables and functions used for tree reconstrcution from incomplete matrices  */

double	**Dist, **Da, *Long, Dmax=0., eps=0.001;
char	**Et, Nom[20], car;
int		N, Ns, NbInc=0, *Tree, *Ch, *Pl, *Clx, *Cly;
int miss=-99;

/* Function prototypes */

void GL_Main(double **DD, int nbObj, int methode, double **distanceArbre, double *RESULTATS, long int *ARETE, double *LONGUEUR, double **W);		
void LecDist();
void EditDist(double **D, int N);
void Chemin(int *T, int t, int u);
void Arete(double *path, int *x, int *y);
int GrowTree(int s, int t, int u, int NbPl);
void MajDist(int s, int NbPl);
void IniTri(double **D, int *ps, int *pt, int *pu);
int PSM(double **D, int N);
void alea (int T, int *aa);
void Ultra1(int *sum, int *sum2, int T, int **B2, double **m, double *max, int *aa);
void Ultra2(int *sum, int *sum2, int T, int **B2, double **m, double *max);
void Additif(int *sum, int *sum2, int T, int **B2, double **m, double *max);		
void compute_criteres11 (double **D,double **DA, double *R, double *LONGUEUR, int nn);
int constructionMisVal (double **D,double **DA,int *X,double **W); 
int parcour211(double **DISS,double **W,double **TM,int *Iternum,long int *ARETE, double *LONGUEUR);
int parcour2MisVal(double **DISS,double **W,double **TM,int *Iternum,long int *ARETE, double *LONGUEUR);
int seqToDistance(char *inputfile,double ***matDistances,int ecrireResultats);
int trex(char*input,char*treeFile,char*statFile,long int **lgintARETE,double **dbLONGUEUR,int ecrireResultat, int unproc);
int bootstrap_function();


/*********************************************************
 * Cette fonction permet de formater le fichier d'entree *
 * au format Phylip, et d'effectuer certaine validation. *
 *********************************************************/

void formatPhylip( char *input){
	int nbSpecies,taille;
	FILE *in,*out;
	char car,car1;
	int newline = 0,i,compt,decompt,retour;   
	char tmp[100];   
	char name[20];
	int nbcar;
	
	get_path(input,tmp);
	strcat(tmp, "tmp");

	in = fopen(input,"r");  
	out = fopen(tmp,"w+"); 
	
	/* Quelques controls d'erreurs*/
	if(m_dataType == 2){
		retour = fscanf(in,"%d %d",&nbSpecies,&taille);
		if(retour<2) { printf("Bad file format !<br>"); exit(8);}
		fprintf(out,"%d %d\n",nbSpecies,taille);
	}
	if(m_dataType == 1){
		retour = fscanf(in,"%d",&nbSpecies);
		if(retour==-1) { printf("Bad file format !<br>"); exit(8);}
		fprintf(out,"%d\n",nbSpecies);
	}
	if(nbSpecies <= 1) { printf("Invalid species number !<br>"); exit(8);}
	
	do{
		car = fgetc(in);     
	}while(car != '\n');
	
	nbcar=11;
	for(i=1;i<=nbSpecies;i++){
		if(nbcar<=10)
			{ printf("Bad file format !<br>"); exit(8);}
		else
			nbcar = 0;
		compt=0;
		
		if(feof(in)) 
			{ printf("Bad file format !<br>"); exit(8);}
		
		do{
			car = fgetc(in); 
			nbcar++;
			if(car==-1) break;
			if(compt < 10){
				if(car == ' ') name[compt] = '_';
				else name[compt] = car;
				compt++;                 
			} 
			else if(compt == 10){
				decompt = 9;
				name[decompt+1] = '\0';
				do{
					if(name[decompt] == '_'){
						name[decompt] = ' ';
					}
				}while(name[--decompt] == '_');
				fprintf(out,"%s ",name);
				fputc(car,out);
				compt++;
			}
			else{
				fputc(car,out);
			}
		}while((car!='\n')&&(car!=EOF));
	}
	while(!feof(in)){
		car = fgetc(in);
		if(!feof(in))
			fputc(car,out);
	} 

	fclose(in);
	fclose(out);
	rename(tmp,input);
}  /* END formatPhylip*/


/*******************
 *La fonction MAIN *
********************/
int main (int argc,char** argv)			
{
	int retour = 0;
	long int *ARETE;
	double *LONGUEUR,**matDist;
	
	/*Initialisation des variables MPI*/
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Get_processor_name(processor_name,&namelen);

	LireOptions(argv,argc);

	//formatPhylip(p_input);

	if((p_optionTR == 1)||(p_optionTR == 2)){ /*inférence d'arbres ou de réticulogramme à partir de matrice*/
		if((m_dataType == 2)&&(p_boot == 2)){ /* sequences sans bootstrap*/
			retour = seqToDistance(p_input,&matDist,1);
			if(retour==0){
				retour = trex(p_matrix,p_output,p_stat,&ARETE,&LONGUEUR,1,0);
			}
		}
		if((m_dataType == 2)&&(p_boot < 2)){ /* sequences avec bootstrap*/
			retour = bootstrap_function();
		}
		if(m_dataType == 1){
			retour = trex(p_input,p_output,p_stat,&ARETE,&LONGUEUR,1,0);
		}
	}
	/* Matrices incompletes */
	else{
		retour = trex(p_input,p_output,p_stat,&ARETE,&LONGUEUR,1,0);
	}

	MPI_Finalize();
	
	return retour;
}/* END main*/


/***********************
 * Recherche du chemin *
 ***********************/
void get_path(char *str,char *path)
{
	register int i;
	
	strcpy(path,str);
	for(i=strlen(path)-1;i>-1;--i) {
		if(str[i]==DIRDELIM) {
			i = -1;
			break;
		}
		if(str[i]=='.') break;
	}
	if(i<0)
		strcat(path,".");
	else
		path[i+1]=EOS;
} /* END get_path*/



/* Récupération de la BIpartition */
void p_RecupererBipartition(long int *ptabARETE, int **pmatBiPartition, int nn)
{
	int i, j, k, l, m, intNbreNoeuds;
	int blnAreteD, blnAreteF;
	long int NoeudDroite, NoeudGauche;
	/* On a n-3 arêtes internes */
	long int *AreteDeb = (long int*)malloc(((nn-3)+1)*sizeof(long int));
	long int *AreteFin = (long int*)malloc(((nn-3)+1)*sizeof(long int));
	long int *ListNInternesD = (long int*)malloc(((nn-2)+1)*sizeof(long int));
	long int *ListNInternesG = (long int*)malloc(((nn-2)+1)*sizeof(long int));
	
	for (i=1; i<=2*nn-3; i++)
		for (j=1; j<=2*nn-3; j++)
			pmatBiPartition[i][j] = 0;
	
	
	/* On sauvegarde les arêtes internes*/
	j = 1;
	for (i=0; i<2*(2*nn-3); i++)
	{	if ((ptabARETE[i] > nn) && (ptabARETE[i+1] > nn)) /* arête interne*/
		{	AreteDeb[j] = ptabARETE[i];
			AreteFin[j] = ptabARETE[i+1];
			j++;
		}
		i++; /* Pour travailler sur 2 cases après 2*/
	}
	
	
	for (i=0; i<2*(2*nn-3); i++)
	{	if (ptabARETE[i] <= nn)  /* une feuille*/
		pmatBiPartition[i/2 + 1][ptabARETE[i]] = 1;
		else if (ptabARETE[i+1] <= nn) /* une feuille*/
			pmatBiPartition[i/2 + 1][ptabARETE[i+1]] = 1;
		else
		{	for (j=1; j<=nn-2; j++)
			{	ListNInternesD[j] = 0;
				ListNInternesG[j] = 0;
			}
			k = 1; l = 1;
			ListNInternesD[k] = ptabARETE[i];
			ListNInternesG[l] = ptabARETE[i+1];
			intNbreNoeuds = (nn-2) - 2; /* nn-2 = Nombre de noeuds internes*/
			while (intNbreNoeuds > 0)
			{	for (j=1; j<=nn-3; j++) /* Parcours de toutes les arêtes internes*/
				{	
					/* Si c'est l'arête pour laquelle on fait cette recherche, 
					   on passe à l'arête suivante*/
					if ((AreteDeb[j] == ptabARETE[i]) &&
						(AreteFin[j] == ptabARETE[i+1])) continue;
					
					/* Récupérer tous les noeuds internes qui sont à droite
					   du noeud droit de l'arête en question*/
					blnAreteD = blnAreteF = 0;
					for (m=1; m<=k; m++)
					{	if (AreteDeb[j] == ListNInternesD[m]) blnAreteD = 1;
						if (AreteFin[j] == ListNInternesD[m]) blnAreteF = 1;
					}
					if (blnAreteD && !blnAreteF)
						{	k++; ListNInternesD[k] = AreteFin[j]; intNbreNoeuds--;	}
					else if (!blnAreteD && blnAreteF)
						{	k++; ListNInternesD[k] = AreteDeb[j]; intNbreNoeuds--;	}
					
					
					/* Récupérer tous les noeuds internes qui sont à gauche
					   du noeud gauche de l'arête en question*/
					blnAreteD = blnAreteF = 0;
					for (m=1; m<=l; m++)
					{	if (AreteDeb[j] == ListNInternesG[m]) blnAreteD = 1;
						if (AreteFin[j] == ListNInternesG[m]) blnAreteF = 1;
					}
					if (!blnAreteD && blnAreteF)
						{	l++; ListNInternesG[l] = AreteDeb[j]; intNbreNoeuds--;	}
					else if (blnAreteD && !blnAreteF)
						{	l++; ListNInternesG[l] = AreteFin[j]; intNbreNoeuds--;	}
				}
			}
			
			/* On va travailler seulement sur les noeud qui se trouvent à droite des noeuds
			   droits des arêtes internes*/
			
			/* Parcours des arêtes depuis le début   */
			
			for (k=0; k<2*(2*nn-3); k++)
			{	NoeudDroite = ptabARETE[k];
				NoeudGauche = ptabARETE[k+1];
				
				/* On cherche les feuilles (<=nn) qui sont liées aux noeuds stoqués
				   dans la liste ListNInternesD, et affecte 1 à la case correspondante
				   dans la matrice des bipartitions*/
				if ((NoeudDroite <= nn) || (NoeudGauche <= nn))
				{
					if (NoeudDroite <= nn)
					{	for (l=1; l<=nn-2; l++)
						if (ListNInternesD[l] == NoeudGauche)
							pmatBiPartition[i/2 + 1][NoeudDroite] = 1;
					}
					else
					{	for (l=1; l<=nn-2; l++)
						if (ListNInternesD[l] == NoeudDroite)
							pmatBiPartition[i/2 + 1][NoeudGauche] = 1;
					}
				}
				k++; /* Pour travailler sur 2 cases après 2*/
			}
		}
		
		i++; /* Pour travailler sur 2 cases après 2*/
	}
	
	free(AreteDeb);
	free(AreteFin);
	free(ListNInternesD);
	free(ListNInternesG);
} /* END p_RecupererBipartition */

/* Comparaison entre les BIpartitions */ 
void p_ComparerBipartition(int **pmatBiPartition, int **pmatBiPartitionTmp, int *ptabComparaison, int nn)
{
	int i, j, k;
	int blnExact, blnExactInvers;
	
	for (i=1; i<=2*nn-3; i++)
	{	for (j=1; j<=2*nn-3; j++)
		{	blnExact = 1; blnExactInvers = 1;
			for (k=1; k<=nn; k++)
			{
				if (blnExact)
					if (pmatBiPartition[i][k] != pmatBiPartitionTmp[j][k]) blnExact = 0;
				if (blnExactInvers)
					if (pmatBiPartition[i][k] != !pmatBiPartitionTmp[j][k]) blnExactInvers = 0;
			}
			if (blnExact || blnExactInvers) 
				{	ptabComparaison[i]++; break; }
		}
	}
} /* END p_RecupererBipartition */

/*******************************************************
//  Algorithme bootstrap                               *
//	1. sequence	-> matrice de distance (D)             *
//  2. matrice de distance -> arbre A(arete,longueur)  *
//  3. bipartition(arete)                              *
//  ........                                           *
********************************************************/
int bootstrap_function(){
	
	/* Calcul pour la séquence initiale*/
	long int *ARETE;
	int *tabComparaison;
	int **TtabComparaison; /* Matrice de comparaison necessaire pour le parallelisme */
	double *LONGUEUR,Pr;
	int intNbreEspeces,intNbreSites,i,r,intBootstrap,intRetour=0;
	char **matSites;	
	double **matDistances;
	char **toto;
	int intDimTab;
	double *m_DI;
	int ** matBiPartition;
	int ** matBiPartitionTmp; 
	char **matNomsEspeces;
	int j,k,p,m = 0;    

	int tour;
	FILE *bs;
	int blnExact, blnExactInvers;
	
	FILE * tmp = fopen(p_input,"r");
	fscanf(tmp,"%d",&intNbreEspeces);
	fscanf(tmp,"%d",&intNbreSites);
	fclose(tmp);
	if(intNbreEspeces<numprocs)
		maxprocs=intNbreEspeces;
	else
		maxprocs = numprocs;
	seqToDistance(p_input,&matDistances,1);
	/* Appel de la fonction trex en mode multiprocesseur */
	trex(p_matrix,"","",&ARETE,&LONGUEUR,0,0);     
	
	matBiPartition = (int **) malloc(((2*intNbreEspeces-3)+2)*sizeof(int*));
	matBiPartitionTmp = (int **) malloc(((2*intNbreEspeces-3)+2)*sizeof(int*));
	for (i=0;i<=2*intNbreEspeces-3+1;i++){	
		matBiPartition[i]=(int *) malloc(((2*intNbreEspeces-3)+2)*sizeof(int));
		matBiPartitionTmp[i]=(int *) malloc(((2*intNbreEspeces-3)+2)*sizeof(int));
	}
	for (i=1; i<=2*intNbreEspeces-3; i++)
		for (j=1; j<=2*intNbreEspeces-3; j++)
			matBiPartitionTmp[i][j] = 0;
	matNomsEspeces=(char **)malloc((intNbreEspeces + 1)*sizeof(char*));
	for (i = 0; i <= intNbreEspeces; i++)
	{	/* 20 correspend à la taille max du nom d'espèce*/
		matNomsEspeces[i] = (char*)malloc((20)*sizeof(char));
	}
	
	if (intNbreSites != 0) /* Matrice de sites*/
		matSites=(char **)malloc((intNbreEspeces +1)*sizeof(char*));	
	matDistances = (double **)malloc((intNbreEspeces + 1)*sizeof(double*));
	for (i = 0; i <= intNbreEspeces; i++)
	{	matSites[i] = (char*)malloc((intNbreSites + 1)*sizeof(char));
		matDistances[i] = (double*)malloc((intNbreEspeces + 1)*sizeof(double));
	}

	p_RecupererBipartition(ARETE,matBiPartition,intNbreEspeces);
	
	strcpy(infilename, p_input);
	infile = fopen(infilename,"r");
	
	fseek(infile, 0, 0);
	get_path(infilename,outfilename);
	strcat(outfilename, "out");
	ibmpc = IBMCRT;
	ansi = ANSICRT;
	if (p_boot == 0)
	{	bootstrap = 1;
		jackknife = 0;
		intBootstrap = 1;
	}
	else
	{	bootstrap = 0;
		jackknife = 1;
		intBootstrap = 2;
	}
	doinput(1, toto);
	
	reps = p_nbRep;

	if (outfilename[0] == '\0')
		{	printf("Error output file.."); return 8; }
	else
	{	
		if(numprocs>1){
			TtabComparaison = (int **) malloc((numprocs+1)*sizeof(int));
			for (r=0;r<=numprocs;r++){
				TtabComparaison[r]=(int *)malloc(((2*intNbreEspeces-3) + 1)*sizeof(int)); 
				for(tour = 1; tour <= 2*intNbreEspeces-3; tour++) 
					TtabComparaison[r][tour] = 0;
			}
		}
		/* Allocation de mémoire du vecteur des pourcentages des arêtes de bases*/
		tabComparaison = (int *)malloc(((2*intNbreEspeces-3) + 1)*sizeof(int)); 
		for (r = 1; r <= 2*intNbreEspeces-3; r++) 
			tabComparaison[r] = 0;
		
		/* Cas de l'exécution sur un seul processeur */ 
		if(numprocs==1){
			for (r = 1; r <= reps; r++) {	
				totota = r;
				Pr=floor((double)(1000 * (r+1)) / (reps + 2));
				
				bootwrite();
				
				fseek(outfile, 0, 0);
				fscanf(outfile,"%d",&intNbreEspeces);
				fscanf(outfile,"%d",&intNbreSites);
				fclose(outfile);
				
				intRetour = seqToDistance(outfilename,&matDistances,1);
				
				if (intRetour != 0)
					return -1;
				/* Appel de la fonction trex en mode uniprocesseur */
				trex(p_matrix,"","",&ARETE,&LONGUEUR,0,1);
				
				p_RecupererBipartition(ARETE, matBiPartitionTmp, intNbreEspeces);
				
				p_ComparerBipartition(matBiPartition, matBiPartitionTmp, tabComparaison, intNbreEspeces);

				if (r < reps)
				{	if ((outfile = fopen(outfilename,"w+")) == 0)                        
					{	printf("Error output file.."); return -1; }
				}
			}
		}
		else{
			for (r = 1; r <= reps; r++) {
				/* Répartition des tâches au niveau de chaque processeur */
				if(r%numprocs==myid){
					totota = r;
					Pr=floor((double)(1000 * (r+1)) / (reps + 2));
					bootwrite();	
					fseek(outfile, 0, 0);
					fscanf(outfile,"%d",&intNbreEspeces);
					fscanf(outfile,"%d",&intNbreSites);
					fclose(outfile);				
					intRetour = seqToDistance(outfilename,&matDistances,1);	
					if (intRetour != 0)
						return -1;
					/* Appel de la fonction trex en mode uniprocesseur */
					trex(p_matrix,"","",&ARETE,&LONGUEUR,0,1);
	
					p_RecupererBipartition(ARETE, matBiPartitionTmp, intNbreEspeces);
	
					p_ComparerBipartition(matBiPartition, matBiPartitionTmp, TtabComparaison[myid], intNbreEspeces);
					
					if (r < reps)
					{	if ((outfile = fopen(outfilename,"w+")) == 0)                        
						{	printf("Error output file.."); return -1; }
					}
				}
			}

			/*Synchronisation du résultat*/
			if(myid>0){
				/* Envoie de matrice de comparaison */
				for (r = 1; r <= 2*intNbreEspeces-3; r++)
					MPI_Send( &TtabComparaison[myid][r], 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
				/* Reception de table de comparaison */
				for (r = 1; r <= 2*intNbreEspeces-3; r++)
					MPI_Recv( &tabComparaison[r], 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				/* Envoie de la table de comparaison au processeur suivant */
				if(myid<numprocs-1){
					for (r = 1; r <= 2*intNbreEspeces-3; r++)
						MPI_Send( &tabComparaison[r], 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );	
				}

			}
			else{
				/* Reception de la matrice de comparaison envoyée par chaque processeur */
				for(tour=1; tour<numprocs;tour++){
					for (r = 1; r <= 2*intNbreEspeces-3; r++)
						MPI_Recv( &TtabComparaison[tour][r], 1, MPI_INT, tour, 0, MPI_COMM_WORLD, &status );
				}
				/* Reconstruction du table de comparaison */
				for (r = 1; r <= 2*intNbreEspeces-3; r++){
					for(tour=0; tour<numprocs;tour++)
						tabComparaison[r] = TtabComparaison[tour][r] + tabComparaison[r];
				}
				/* Envoie du tableau de comparaison au processeur suivant */
				for (r = 1; r <= 2*intNbreEspeces-3; r++)
					MPI_Send( &tabComparaison[r], 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );	
			}
		}

	}
	/* Calcul des pourcentages des arêtes internes*/
	for (i=1;i<=2*intNbreEspeces-3;i++)
		tabComparaison[i] = tabComparaison[i] * 100 / reps;
	
	seqToDistance(p_input,&matDistances,1);


	trex(p_matrix,"","",&ARETE,&LONGUEUR,0,0);
	bs = fopen(p_outboot,"w+");
	for (i=1;i<=2*intNbreEspeces-3;i++){
		fprintf(bs,"%d %d %d\n",ARETE[2*i-1],ARETE[2*i-2],tabComparaison[i]);
	}
	fclose(bs);

	trex(p_matrix,p_output,p_stat,&ARETE,&LONGUEUR,1,0);

	return 0;
} /* END bootstrap_function */


/* Tree reconstruction from sequences*/
int seqToDistance(char *inputfile,double ***matDistances_,int ecrireResultats){
	int intNbreEspeces;
	int intNbreSites;
	int i,intRetour;
	char **matNomsEspeces, **matSites;
	double **matDistances;
	
	int m_intIdMethod	= m_sequenceMethod;
	double m_dblPenalty = m_gapPenalty; 
	double 	m_dblParamA = m_paramA;
	
	FILE *fileSource;
	
	/*CalculDistances* instCalcDist = new CalculDistances;*/
	fileSource = fopen(inputfile,"r");
	
	fscanf(fileSource,"%d",&intNbreEspeces);
	fscanf(fileSource,"%d",&intNbreSites);
	
	matNomsEspeces=(char **)malloc((intNbreEspeces + 1)*sizeof(char*));
	for (i = 0; i <= intNbreEspeces; i++)
	{	/* 20 correspend à la taille max du nom d'espèce*/
		matNomsEspeces[i] = (char*)malloc((20)*sizeof(char));
	}
	
	if (intNbreSites != 0) /*Matrice de sites*/
		matSites=(char **)malloc((intNbreEspeces +1)*sizeof(char*));	
	matDistances = (double **)malloc((intNbreEspeces + 1)*sizeof(double*));
	for (i = 0; i <= intNbreEspeces; i++)
	{	matSites[i] = (char*)malloc((intNbreSites + 1)*sizeof(char));
		matDistances[i] = (double*)malloc((intNbreEspeces + 1)*sizeof(double));
	}
		
	intRetour = f_CalculerDistances(fileSource, p_matrix,matNomsEspeces, matDistances, matSites,intNbreEspeces, intNbreSites,m_intIdMethod,m_dblPenalty, 	m_PEMV,m_dblParamA,ecrireResultats,0);
	fclose(fileSource);
	
	if (intRetour != 0)
		return -1;
	
	(*matDistances_) = matDistances;
	
	return 0;
} /* END seqToDistance */

void SAVEASNewick(char *string, double *LONGUEUR, long int *ARETE, char **Et, int nn) 
{
	int n,root,a;
	int Ns;
	int i, j, sd, sf, *Suc, *Fre, *Tree, *degre, *Mark;
	double *Long; 
	/*if (nn!=NULL)*/
	if (nn !=0) n = nn;
	Ns=2*n-3;
	
	Suc =(int*) malloc((2*n) * sizeof(int));
	Fre =(int*) malloc((2*n) * sizeof(int));
	degre =(int*) malloc((2*n) * sizeof(int));
	Long = (double*) malloc((2*n) * sizeof(double));	
	Tree = (int*) malloc((2*n) * sizeof(int));
	Mark =(int*) malloc((2*n) * sizeof(int));
	
	if ((degre==NULL)||(Mark==NULL)||(string==NULL)||(Suc==NULL)||(Fre==NULL)||(Long==NULL)||(Tree==NULL)||(ARETE==NULL)||(LONGUEUR==NULL))	
		{ printf("Tree is too large to be saved"); return;} 
	
	for (i=1;i<=2*n-3;i++)
	{ 
		
		if (i<=n) degre[i]=1;
		else degre[i]=3;
	} degre[2*n-2]=3;
	
	root=Ns+1;
	for (;;)
	{         
		a=0; a++;
		for (j=1;j<=2*n-2;j++)
			Mark[j]=0;
		
		for (i=1;i<=2*n-3;i++)
		{ 	  									
			if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
				
			}
			else if ((degre[ARETE[2*i-1]]==1)&&(degre[ARETE[2*i-2]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-1]]=ARETE[2*i-2]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-1]]=LONGUEUR[i-1];
				
			}
			else if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]==1)&&(Mark[ARETE[2*i-2]]==0)&&(Mark[ARETE[2*i-1]]==0))
			{
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; root=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; a=-1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
			}
			if (a==-1) break;
		}
		if (a==-1) break;
	}
	
	
	/*  On decale et on complete la structure d'arbre avec Successeurs et Freres  */
	for (i=Ns+1;i>0;i--)
	{ 	Fre[i]=0; Suc[i]=0;
	}	Tree[root]=0;/*Tree[Ns+1]=0;*/
	
	for (i=1;i<=Ns+1/*Ns*/;i++)
	{	
		if (i!=root) 
		{
			sd=i; sf=Tree[i];
			if (Suc[sf]==0) Suc[sf]=sd;
			else {	
				sf=Suc[sf];
				while (Fre[sf]>0) sf=Fre[sf];
				Fre[sf]=sd;
			}		 
		}
	}
	
	
	/* On compose la chaine parenthesee */
	string[0]=0; i=root;/*i=Ns+1;*/
	for (;;)
	{	if (Suc[i]>0)
		{	sprintf(string,"%s(",string);
		Suc[i]=-Suc[i]; i=-Suc[i]; }
		else if (Fre[i]!=0)
		{	if (Suc[i]==0) sprintf(string,"%s%s:%.4f,",string,Et[i-1],Long[i]);
			else sprintf(string,"%s:%.4f,",string,Long[i]);
		i=Fre[i]; }
		else if (Tree[i]!=0)
		{	if (Suc[i]==0) sprintf(string,"%s%s:%.4f)",string,Et[i-1],Long[i]);
		else sprintf(string,"%s:%.4f)",string,Long[i]); i=Tree[i]; }
		else break;
	}	
	strcat(string,";");
	free(Suc); free(Fre); free(Tree); free(Long); free(degre); free(Mark);	
} /* END SAVEASNewick */



/*Fonction principale du logiciel trex*/
/* unproc = 0(non) et 1(oui)*/
	
int trex(char *inputFile,char *treeFile, char *statFile,long int **ARETE_,double **LONGUEUR_,int ecrireResultats, int unproc){ 
	double Puissance;
	int l,i,j,Iternum,ReticulationsNumber=0, OptionFunction=0, optionMiss=0,method=0,optionTR=0;
	double *LONGUEUR,power;
	long int *ARETE;
	double **DI,**DA,**DA1,**W,*R,*RAJ; 
	int Option,Option1;
	double *CRITERION1, *CRITERION2;
	double RESULTATS[4];
	char mwFile[50];
	FILE *mod;
	FILE *data;
	FILE *newickFile;
	
	char /*MW[20],*/ symbol, PrintTreeMetric='Y';
	char **Names;
	double c;
	int k;
	char car;
	char *stringNewick = (char*)malloc(1000);
	char chCaractere;
	
	if ((data=fopen(inputFile,"r"))==0) { printf("File %s was not found \n", p_input); exit(8); } 
	
	/* Reading the dissimilarity data file */
	fscanf(data,"%d",&n);
	if ((n>10000)||(n<1)) { printf("*Incorrect Phylip format or too many taxa in the distance matrix n = %d\n\n", n); exit(5); } 
	if(n<numprocs)
		maxprocs=n;
	else
		maxprocs = numprocs;
	/*Merory allocation*/
	DI=(double **) malloc((n+1)*sizeof(double*));
	DA=(double **) malloc((n+1)*sizeof(double*));
	DA1=(double **) malloc((n+1)*sizeof(double*));
	W=(double **) malloc((n+1)*sizeof(double*));
	R=(double *) malloc(4*sizeof(double));
	RAJ=(double *) malloc(4*sizeof(double));
	ARETE=(long int *) malloc((4*n-1)*sizeof(long int)); 
	LONGUEUR=(double *)malloc((2*n-1)*sizeof(double));
	Names=(char **)malloc((n+1)*sizeof(char*));
	
	for (i=0;i<=n;i++)
	{
		DI[i]=(double*)malloc((n+1)*sizeof(double));
		DA[i]=(double*)malloc((n+1)*sizeof(double));
		DA1[i]=(double*)malloc((n+1)*sizeof(double));
		W[i]=(double*)malloc((n+1)*sizeof(double));
		Names[i]=(char *)malloc((50)*sizeof(char));
		
		if ((DI[i]==NULL)||(DA[i]==NULL)||(DA1[i]==NULL)||(W[i]==NULL)||(Names[i]==NULL))
		{
			printf("*Data matrix is too large\n "); 
			exit(5);
		}
	}

	/* Reading the dissimilarity data file */
	for (i=1; i<=n; i++)
	{
		while (((symbol=getc(data))==' ')||(symbol=='\n')||(symbol=='\t'));
		fseek(data,-1,1);
		k=0;
        while (((symbol=getc(data))!=' ')&&(symbol!='\t'))
        {
            if (k<30) Names[i-1][k] = symbol;
            k++;
        }
        Names[i-1][k] = '\0';

        for (j=1; j<=n; j++)
            fscanf(data,"%lf",&DI[i][j]);
	}

	DI[0][0]=0;
	fclose (data);
	optionTR=p_optionTR;
	/*Reconstruction method:
	1. Infer an additive (phylogenetic) tree from D or Infer an additive tree from sequences
	2. Infer a reticulogram (reticulated network) from D
	3. Infer an additive tree from an incomplete matix D*/

	
	/*Tree reconstruction from incomplete matrices*/
	optionMiss=p_method;
	if (optionTR==3){
		optionMiss=p_method;
		/*  Reconstruction method for an incomplete dissimilarity matrix D\n(missing entries in D should be indicated by -99)
			1. Triangles method - Guenoche, Leclerc (2001");
			2. Ultrametric procedure + MW - De Soete (1984));
			3. Additive procedure + MW - Landry et al. (1996));
			4. MW-modified - Makarenkov (2001);
			5. MW* - Makarenkov, Lapointe (2004)*/
		GL_Main(DI, n, optionMiss, DA, RESULTATS, ARETE, LONGUEUR, W);	
	} 
	/* Tree reconstruction from comlete matrices and reticulogram reconstruction*/
	else {
		
		method=p_method;		
	    /* Tree reconstruction menu (5 methods available) :
		   1. ADDTREE - Sattath and Tversky (1977)
		   2. Neighbor Joining - Saitou and Nei (1987)
		   3. Unweighted Neighbor Joining - Gascuel (1997)
		   4. Circular order reconstruction - Makarenkov, Leclerc (1997)
		   5. Weighted least-squares method MW - Makarenkov, Leclerc (1999)*/
		
		if (method==5) /* options for the MW method*/
		{
			if (n>3)
			{
				Option = p_option;
			    /*The MW optimization options :
				  1. Global MW optimization
				  2. Local MW optimization*/
				DI[0][0]=Option;   
			}
			else 
				Option=1;

			Option1 = p_option1;
		    /*Choice of the weight matrix W :
				1. Of the form W=1/D^p
				2. In a file */
			
			if (Option1==1) 
			{
				Puissance = p_p;
				for (j=1;j<=n;j++)
				{
					for (i=j+1;i<=n;i++)
					{       
						if (DI[j][i]>0)
							W[j][i]=pow(1/DI[j][i],Puissance); 
						else if (DI[j][i]<=0)    
							W[j][i]=0.000001;                      
						W[i][j]=W[j][i];
					}
					W[j][j]=1;
				}  
			}
			
			else if (Option1==2) /* reading W from a file*/
			{   
				FILE *weights;   
				if ((weights=fopen(mwFile,"r"))==0) { printf("*File %s was not found ", mwFile); 
				exit(8); } 
				for (j=2; j<=n; j++)
				{
					for (i=1; i<=j-1; i++)
					{
						fscanf(weights,"%lf",&c);
						W[i][j]=c;
						W[j][i]=c;
					}
					W[j][j]=0;
				}
				fclose (weights); 	
			}
		}
		
		if (method!=5) /* Setting up W matrix if MW is not selected*/
		{
			for (i=1;i<=n;i++)
			{
				for (j=i+1;j<=n;j++)
				{
					W[i][j]=1;
					W[j][i]=W[i][j];
				}
				W[i][i]=1;
			}
		}      
		W[0][0]=1;
		
		if (n==2) /* the trivial case - of two objects*/
		{
			DA[1][2]=DI[1][2];
			DA[2][1]=DI[1][2];
			DA[1][1]=0.0;
			DA[2][2]=0.0;
		} 
		
		else 
		{    
			if (n<=10) l=30;
			else l=3*n;
			Iternum=10;
			if (optionTR==1)
				Scaling(DI,&power);
			
			/* Blocage pour attendre les autres processeurs */
			if(unproc==0)
				MPI_Barrier(MPI_COMM_WORLD);

			if  (method==1)  /*method ADDTREE*/
				ADDTREE(DI,DA1,unproc);  

			else if  (method==2) /*method NJ*/
				NJ(DI,DA1,unproc);

			else if (method==4)  /*method Circular Orders*/
				parcour1(DI,W,DA,&Iternum,ARETE,LONGUEUR,unproc);
			
			else if  (method==3)  /*method UNJ*/
				UNJ(DI,DA1,unproc);  
		      
			else if (method==5) /*method MW*/
			{   		   		   		
				if (Option==1)
					DI[0][0]=1;
				else 	
					DI[0][0]=2;
				parcour2(DI,W,DA,&Iternum,ARETE,LONGUEUR,unproc);
			}
			else if(method==6)/*method bioNJ*/
				BioNJ_main (DI,DA);
				
			if (method<4 || method==6) /*edge lengths polishing for ADDTREE, NJ, and UNJ*/
				approx_arb2 (DI,DA1,DA,W,&l,ARETE,LONGUEUR);     
		}  
		
		/* Reticulogram reconstruction option */
		if (optionTR==2) 
		{ 
			int Iternumber=(2*n-2)*(2*n-2-1)/2-2*n+3;  
			
			/*Selection of the stopping rule*/
			OptionFunction = p_optionFunction;
			
			/*	Stopping rule for addition of new edges to the reticulogram:
				1. When the criterion Q1 is minimized.
				2. When the criterion Q2 is minimized.
				3. Add a fixed number of edges K to the reticulogram.*/
			
		    /*if (OptionFunction==3)
			{ 
				Iternumber = p_iternumber;
				while ((Iternumber<0)||(Iternumber>=(2*n-2)*(2*n-2-1)/2-2*n+3)){
					printf("Number of edges must be between 0 and %d\n\n",((2*n-2)*(2*n-2-1)/2-2*n+3)-1);  
					exit(8);
					printf("\n\n");
				}
			}  */
			CRITERION2 = (double*)malloc((2*n-3 + Iternumber+1)*sizeof(double));
			CRITERION1 = (double*)malloc((2*n-3 + Iternumber+1)*sizeof(double)); 
			
			free(ARETE); free(LONGUEUR);
			ARETE=(long int *) malloc(((2*n-2)*(2*n-3)+1)*sizeof(long int)); 
			LONGUEUR=(double *)malloc(((2*n-2)*(2*n-3)/2+1)*sizeof(double));	
		/*	if(unproc==0)
				MPI_Barrier(MPI_COMM_WORLD);*/
			RETICULATIONS (n, DI, DA, ARETE, LONGUEUR, OptionFunction, W, Iternumber, &ReticulationsNumber, CRITERION1, CRITERION2, unproc);		
		}
		else{
			CRITERION2 = (double*)malloc(sizeof(double)); 
			CRITERION1 = (double*)malloc(sizeof(double)); 
		}
	}

	/* Écriture des résultats pour le premier processeur seulement */
	if(ecrireResultats && myid==0){
		printf("\n %s", stringNewick);
		SAVEASNewick(stringNewick,LONGUEUR, ARETE, Names,n) ;

		newickFile = fopen(p_newick,"w+");
		fprintf(newickFile,"%s",stringNewick);
		fclose(newickFile);
	
		TreeFile = fopen(treeFile,"w+");
		if (TreeFile==NULL) {printf("\nFichier %s inexistant xxxxxxxxxxxxxxxx",treeFile); exit(8);}
		
	
		mod = fopen(p_modele,"r");
		if (mod==NULL) {printf("\nFichier %s inexistant",p_modele); exit(8);}
		
		car = fgetc(mod);
		do{		
			fputc(car,TreeFile);
			car = fgetc(mod);
		}while (car!=EOF);
		fclose(mod);
		fprintf(TreeFile,"drawRetic = 1\n");
		fprintf(TreeFile,"n = %d\n",n);
		
	
		Output1=fopen(statFile,"w+");
		if(Output1 == NULL) exit(9);
		
		if ((method==5)&&(Option1==1))  
 			fprintf(Output1,"Weight matix is of the form W=1/D^%f \n",Puissance);
		


		if (optionTR!=3)
		{
			if (method==1) fprintf(Output1,"\nTree reconstruction method - ADDTREE\n"); 
			else if (method==2) fprintf(Output1,"\nTree reconstruction method - Neighbor Joining\n"); 
			else if (method==3) fprintf(Output1,"\nTree reconstruction method - Unweighted Neighbor Joining\n"); 
			else if (method==4) fprintf(Output1,"\nTree reconstruction method - Circular order reconstruction\n"); 
			else if (method==5) 
 			{fprintf(Output1,"\nTree reconstruction method - Weighted least-squares method MW"); 
				if (Option==1) fprintf(Output1," (global optimization)\n");
				else fprintf(Output1," (local optimization)\n");
   			}
			else if (method==6) fprintf(Output1,"\nTree reconstruction method - BioNJ\n"); 
		}
		else 
		{
			fprintf(Output1,"\nTree reconstruction from incomplete dissimilarity matrix D\n\n"); 
			if (optionMiss==1) fprintf(Output1,"\nTree reconstruction method - Triangles method\n"); 
			else if (optionMiss==2) fprintf(Output1,"\nTree reconstruction method - Ultrametric procedure + MW\n"); 
			else if (optionMiss==3) fprintf(Output1,"\nTree reconstruction method - Additive procedure + MW\n"); 
			else if (optionMiss==4) fprintf(Output1,"\nTree reconstruction method - MW-modified\n"); 
			else if (optionMiss==5) fprintf(Output1,"\nTree reconstruction method - MW*\n"); 
		}
		
		if (power!=0.0) W[0][0]=2; 
		else if (PrintTreeMetric=='Y') W[0][0]=0.1;
		else W[0][0]=0;
		
		if (optionTR==3) { RAJ[1]=RESULTATS[0]; RAJ[2]=RESULTATS[1]; RAJ[3]=RESULTATS[2]; RAJ[0]=RESULTATS[3]; }
		compute_criteres(DI,DA,RAJ,W,optionTR,Names);
		
		if ((power!=0.0)&&(n>2))
		{
			if (optionTR==1) Scaling1(DI,DA,RAJ,LONGUEUR,&power);
			if (PrintTreeMetric=='Y') W[0][0]=0.1;
			else W[0][0]=0;
			compute_criteres(DI,DA,RAJ,W,optionTR,Names);
		}
		
		if (n==2) { ARETE[0]=1;ARETE[1]=2;LONGUEUR[0]=DA[1][2]; }
		if (optionTR!=2) ReticulationsNumber=2*n-3;
		PrintEdges(ARETE,LONGUEUR,ReticulationsNumber,optionTR,CRITERION1,CRITERION2,OptionFunction);
		
		fclose (Output1); 
		fclose (TreeFile);
		printf("\n\n");  
	}

	(*ARETE_) = ARETE;
	(*LONGUEUR_) = LONGUEUR;

} /* END trex */

// Method MW

void construction(double **D,double **DA,int *X,double **W)
{
	int i,j,i1,p,P,xi,a,*Y,NV,NR,PP,PN,option,*ARETE;
	double *A,am,c,am1,c1,k,k1,k2,k5,k6,k7,ki,**G1,S1,S2,
	S3,M,MR,M1,DIS,DIS1,Su,Sv,*L,**L1,**L2,
	**L3,*DIST,*DIST1;
	double Dummy1;
	
	
	Y=(int *) malloc((n+1)*sizeof(int));
	A=(double *) malloc((8*n+8)*sizeof(double));
	G1=(double **) malloc((2*n-2)*sizeof(double*));
	
	L=(double *) malloc((n+1)*sizeof(double));
	L1=(double **) malloc((2*n-2)*sizeof(double*));
	L2=(double **) malloc((2*n-2)*sizeof(double*));
	L3=(double **) malloc((2*n-2)*sizeof(double*));
	DIST=(double *) malloc(3*sizeof(double));
	DIST1=(double *) malloc(3*sizeof(double));
	ARETE=(int *) malloc(3*sizeof(int));
	
	
	
	for (i=0;i<=2*n-3;i++)
	{
		G1[i]=(double*)malloc((n+1)*sizeof(double));
		L1[i]=(double*)malloc((n+1)*sizeof(double));
		L2[i]=(double*)malloc((n+1)*sizeof(double));
		L3[i]=(double*)malloc((n+1)*sizeof(double)); 
		if ((G1[i]==NULL)||(L1[i]==NULL)||(L2[i]==NULL)||(L3[i]==NULL))
		{printf("probleme de memoire...");
			exit(5);
		}
	}      
	
	Dummy1=100;
	if (D[X[1]][X[2]]<(n-1)*0.0002){ Dummy1=D[X[1]][X[2]]; D[X[1]][X[2]]=(n-1)*0.0002;}
	
	option=X[0];
	
	MR=0;
	A[1]=X[1];
	A[2]=X[2];
	A[3]=0;
	
	A[4]=D[X[1]][X[2]]; 
	
	for (i=1;i<=n;i++)
		Y[i]=1;
	Y[X[1]]=0;
	Y[X[2]]=0;
	P=1;
	DA[X[1]][X[1]]=0;
	DA[X[2]][X[2]]=0;
	DA[X[1]][X[2]]=A[4];
	DA[X[2]][X[1]]=A[4];
	
	xi=X[2];
	for (j=1;j<=n;j++)
	{
		if ((j!=X[1])&&(j!=X[2]))
		{
			k1=DA[X[1]][xi]-D[xi][j];
			k2=D[X[1]][j];
			L[j]=W[X[1]][j]+W[X[2]][j];
			L1[1][j]=2*(k1*W[X[2]][j]-k2*W[X[1]][j]);
			L2[1][j]=2*(-k1*W[X[2]][j]-k2*W[X[1]][j]);
			L3[1][j]=2*(W[X[1]][j]-W[X[2]][j]);
			G1[1][j]=k1*k1*W[X[2]][j]+k2*k2*W[X[1]][j];
		}
	}
	
	for (i=2;i<=n-1;i++)
	{    a=0;
		for (p=1;p<=P;p++)
			if ((option!=2)||(i>3)) 
		{
			
			for (p=1;p<=P;p++)
			{
				for (j=1;j<=n;j++)
				{
					if (Y[j]==1)
						
					{
						xi=floor1(A[4*p-2]);
						k1=DA[X[1]][xi]-D[xi][j];
						k2=D[X[1]][j];
						Su=A[4*p-1];
						Sv=Su+A[4*p];
						
						
						k6=L1[p][j];
						k5=L2[p][j];
						k7=L3[p][j];
						k=L[j];
						M=G1[p][j]+MR+1.0;
						
						if ((2*k*Sv+k5<=0.00001)&&(k6+k7*Sv>=-0.00001))
						{ M=k*Sv*Sv+k5*Sv;
							am=Sv;
							c=0;
						}
						if (((c1=-(k6+k7*Sv)/(2*k))>=-0.00001)&&(-k5-2*k*Sv-k7*c1>=-0.00001)&&((M1=k*(Sv*Sv+c1*c1)+k7*Sv*c1+k5*Sv+c1*k6)<M))
						{ M=M1;
							am=Sv;
							c=c1;
						}
						if (((am1=-k5/(2*k))>=Su-0.00001)&&(-k6-k7*am1<=0.00001)&&((M1=k*am1*am1+k5*am1)<M)&&(am1<=Sv+0.00001))
						{ M=M1;
							c=0;
							am=am1;
						}
						if (((4*k*k-k7*k7)!=0)&&((am1=(-2*k*k5+k7*k6)/(4*k*k-k7*k7))<=Sv+0.00001)&&((c1=-(am1*k7+k6)/(2*k))>=-0.00001)&&((M1=k*(am1*am1+c1*c1)+k7*am1*c1+k5*am1+k6*c1)<M)&&(am1>=Su-0.00001))
						{
							M=M1;
							c=c1;
							am=am1;
						}
						
						if (((4*k*k-k7*k7)==0)&&(k5==k6)&&((-k5/(2*k))>=Su-0.00001)&&((c1=-k5/(2*k)-Su)>=-0.00001)&&((M1=k*(Su+c1)*(Su+c1)+k5*(Su+c1))<M))
						{
							M=M1;
							c=c1;
							am=Su;
						}
						
						if ((2*k==k7)&&(k5==k6)&&(k5>=-0.00001)&&((M1=k*Su*Su+k5*Su)<M))
						{
							M=M1;
							c=0;
							am=Su;
						}
						
						if ((2*k*Su+k5>=-0.00001)&&(-k6-k7*Su<=0.00001)&&((M1=k*Su*Su+k5*Su)<M))
						{ c=0;
							am=Su;
							M=M1;
						}
						
						if (((c1=-(Su*k7+k6)/(2*k))>=-0.00001)&&((k5+2*k*Su+c1*k7)>=-0.00001)&&((M1=k*Su*Su+k5*Su+c1*c1*k+c1*(Su*k7+k6))<M))
						{ 
							am=Su;
							M=M1;
							c=c1;
						}
						
						if ((a==0)||(MR>M+G1[p][j]))
						{
							a=1;
							MR=M+G1[p][j];
							DIS=am;
							DIS1=c;
							
							if (c<=(n-i+1.0)*0.00002) DIS1=(n-i+1.0)*0.00002;
							if (fabs(Sv-Su)<=2*(n-i+1)*0.00002) DIS=Sv-(Sv-Su)/(n-i+1.0);
							else 
							{   
								if (fabs(am-Su)<=(n-i+1.0)*0.00002) DIS=Su+(n-i+1.0)*0.00002;
								if (fabs(am-Sv)<=(n-i+1.0)*0.00002) DIS=Sv-(n-i+1.0)*0.00002;
							} 
							
							
							NV=j;
							NR=p;
						}
					}
				}
			}
			
		}
		
		if ((option==2)&&(i<=3))
		{
			NV=X[i+1];
			NR=ARETE[i-1];
			DIS=DIST[i-1];
			DIS1=DIST1[i-1];
		}  
		
		
		/*           Arrays  Updating     */
		p=NR;
		PP=P;
		if (fabs(DIS-A[4*p-1])<=0.00001)
			DIS=A[4*p-1];
		if (fabs(DIS-A[4*p-1]-A[4*p])<=0.00001)
		{
			DIS=A[4*p-1]+A[4*p];
		}
		X[i+1]=NV;
		xi=floor1(A[4*p-2]);
		Y[NV]=0;
		if ((DIS-A[4*p-1]>0)&&(DIS-A[4*p-1]-A[4*p]<0))
		{
			c=A[4*p];
			A[4*p]=DIS-A[4*p-1];
			P=P+1;
			A[4*P-3]=P;
			A[4*P-2]=A[4*p-2];
			A[4*P-1]=DIS;
			A[4*P]=A[4*p-1]+c-DIS;
		}
		if (DIS1>0)
		{
			P=P+1;
			A[4*P-3]=P;
			A[4*P-2]=X[i+1];
			A[4*P-1]=DIS;
			A[4*P]=DIS1;
		}
		DA[X[1]][X[i+1]]=DIS+DIS1;
		DA[X[i+1]][X[1]]=DIS+DIS1;
		
		
		/* The induction formula */
		
		for (j=2;j<=i;j++)
		{
			
			if (((DA[X[1]][X[j]]+DA[X[1]][xi]-DA[xi][X[j]])/2)<=DIS)
				DA[X[j]][X[i+1]]=DIS+DIS1+DA[X[j]][xi]-DA[X[1]][xi];
			
			else
				DA[X[j]][X[i+1]]=DIS1-DIS+DA[X[1]][X[j]];
			
			
			DA[X[i+1]][X[j]]=DA[X[j]][X[i+1]];
			DA[X[i+1]][X[i+1]]=0;
			
		}
		
		if ((P==PP+2)||((P==PP+1)&&(DIS1==0)))
			PN=PP+1;
		else
			PN=PP;
		if ((DIS-A[4*p-1]>0)&&(DIS-A[4*p-1]-c<0))
		{
			for(j=1;j<=n;j++)
			{
				if (Y[j]==1)
				{
					L1[PP+1][j]=L1[p][j];
					L2[PP+1][j]=L2[p][j];
					L3[PP+1][j]=L3[p][j];
					G1[PP+1][j]=G1[p][j];
				}
			}
		}
		
		/* BLOCK 1 */
		for (p=1;p<=PN;p++)
		{
			xi=floor1(A[4*p-2]);
			for (j=1;j<=n;j++)
			{
				if (Y[j]==1)
				{ 
					if (((DA[X[1]][xi]+DA[X[1]][X[i+1]]-DA[X[i+1]][xi])/2)<=A[4*p-1]+0.00001)
					{
						ki=D[X[i+1]][j]+DA[X[1]][xi]-DA[X[i+1]][xi];
						L1[p][j]=L1[p][j]-2*ki*W[X[i+1]][j];
						L2[p][j]=L2[p][j]-2*ki*W[X[i+1]][j];
						L3[p][j]=L3[p][j]+2*W[X[i+1]][j];
					}
					else
					{
						ki=D[X[i+1]][j]-DA[X[i+1]][X[1]];
						L1[p][j]=L1[p][j]-2*ki*W[X[i+1]][j];
						L2[p][j]=L2[p][j]+2*ki*W[X[i+1]][j];
						L3[p][j]=L3[p][j]-2*W[X[i+1]][j];
					}
					if (p==1)
						L[j]=L[j]+W[X[i+1]][j];
					G1[p][j]=G1[p][j]+ki*ki*W[X[i+1]][j];
				}
			}
		}
		
		/* BLOCK 2 */
		if (DIS1>0)
		{
			for (j=1;j<=n;j++)
			{
				if (Y[j]==1)
				{
					k1=DA[X[1]][X[i+1]]-D[X[i+1]][j];
					k2=D[X[1]][j];
					S1=W[X[1]][j]+W[X[i+1]][j];
					S2=k2*W[X[1]][j];
					S3=k1*k1*W[X[i+1]][j]+k2*k2*W[X[1]][j];
					for (i1=2;i1<=i;i1++)
					{
						ki=D[X[i1]][j]+DA[X[1]][X[i+1]]-DA[X[i+1]][X[i1]];
						S1=S1+W[X[i1]][j];
						S2=S2+ki*W[X[i1]][j];
						S3=S3+ki*ki*W[X[i1]][j];
					}
					L[j]=S1;
					L1[P][j]=2*(k1*W[X[i+1]][j]-S2);
					L2[P][j]=2*(-k1*W[X[i+1]][j]-S2);
					L3[P][j]=2*(S1-2*W[X[i+1]][j]);
					G1[P][j]=S3;
				}
			}
		}
	}
	
	if (Dummy1<(n-1)*0.0002) D[X[1]][X[2]]=Dummy1;
	
	free(L);
	free(Y);
	free(A);
	free(DIST);
	free(DIST1);
	free(ARETE);
	
	for (i=0;i<=2*n-3;i++)
	{ 
		free(G1[i]);
		free(L1[i]);
		free(L2[i]);
		free(L3[i]);
	}    
	free(G1);
	free(L1);
	free(L2);
	free(L3); 
} /* END construction */

/* Extacting a Circular order X of D*/

void odp(double **D,int *X,int *i1,int *j1)
{
	double S1,S;
	int i,j,k,a,*Y1;

	Y1=(int *) malloc((n+1)*sizeof(int));

	for(i=1;i<=n;i++)
		Y1[i]=1;

	X[1]=*i1;
	X[n]=*j1;
	if (n==2) return;
	Y1[*i1]=0;
	Y1[*j1]=0;
	for(i=0;i<=n-3;i++)
	{ a=2;
	S=0;
	for(j=1;j<=n;j++)
	{ if (Y1[j]>0)
	{
		S1= D[X[n-i]][X[1]]-D[j][X[1]]+D[X[n-i]][j];
		if ((a==2)||(S1<=S))
		{
			S=S1;
			a=1;
			X[n-i-1]=j;
			k=j;
		}
	}
	}
	Y1[k]=0;
	}
	free(Y1);
} /* END odp */

/* Neighbor Joining method */

void NJ(double **D1,double **DA, int unproc)
{
	double **D,*T1,*S,*LP,Som,Smin,Sij,L,Lii,Ljj,l1,l2,l3;
	int *T,i,j,ii,jj,n1;
	double *Tsmin; /* Tableau des Smin necessaire à la parallelisation */
	int *Tii, *Tjj; /* Tableau des ii et jj necessaire à la parallelisation */

	/* Allocation de la mémoire */
	D=(double **) malloc((n+1)*sizeof(double*));
	T1=(double *) malloc((n+1)*sizeof(double));
	S=(double *) malloc((n+1)*sizeof(double));
	LP=(double *) malloc((n+1)*sizeof(double));
	T=(int *) malloc((n+1)*sizeof(int));
	
	/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1 || unproc==1){
		Tii=(int *) malloc((numprocs+1)*sizeof(int));
		Tjj=(int *) malloc((numprocs+1)*sizeof(int));
		Tsmin=(double *) malloc((numprocs+1)*sizeof(double));
	}
	Smin=(double)LONG_MAX;

	/* Allocation de la mémoire */
	for (i=0;i<=n;i++){
		D[i]=(double*)malloc((n+1)*sizeof(double));
		if (D[i]==NULL){
			printf("probleme de memoire");
			exit(5);
		}
	}
	L=0;Som=0;

	/*Initialisation des valeurs*/
	for (i=1;i<=n;i++){
		S[i]=0; LP[i]=0;
		for (j=1;j<=n;j++){
			D[i][j]=D1[i][j];
			S[i]=S[i]+D[i][j];
		}
		Som=Som+S[i]/2;
		T[i]=i;
		T1[i]=0;
	}

	/* Main procedure */
	for (n1=n;n1>3;n1--){
		/* Research of the best pair of objects (i,j) for clustering*/
		Smin=2*Som;
		/* Cas d'un processeur */
		if(numprocs==1 || unproc==1){
			for (i=1;i<=n1-1;i++){
				for (j=i+1;j<=n1;j++){
					Sij=2*Som-S[i]-S[j]+D[i][j]*(n1-2);
					if (Sij<Smin){
						Smin=Sij;ii=i;jj=j;
					}
				}
			}
		}
		/* Cas de plusieurs processeurs */
		else{
			for (i=1;i<=n1-1;i++){
				if(i%numprocs==myid){
					for (j=i+1;j<=n1;j++){
						Sij=2*Som-S[i]-S[j]+D[i][j]*(n1-2);
						if (Sij<Smin){
							Smin = Sij;ii=i; jj=j;Tsmin[myid] = Sij;Tii[myid]=i; Tjj[myid]=j;
						}
					}
				}
			}	
			/* Syncrhonisation des résultats */
			if(myid>0){	
				MPI_Send( &Smin, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
				MPI_Send( &ii, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
				MPI_Send( &jj, 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
				MPI_Recv( &ii, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &jj, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				if(myid<numprocs-1){
					MPI_Send( &ii, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
					MPI_Send( &jj, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				}
			}
			else{
				for(i=1; i<numprocs;i++){
					MPI_Recv( &Tsmin[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Tii[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Tjj[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
				}
				for(i=1; i<numprocs;i++){
					if( Tsmin[i]==Smin && ( ii>Tii[i] || (ii==Tii[i]&& jj>Tjj[i]) )){
						ii = Tii[i];jj = Tjj[i];Smin = Tsmin[i];
					}
					if(Tsmin[i]<Smin){
						ii = Tii[i];jj = Tjj[i];Smin = Tsmin[i];
					}	
				}
				MPI_Send( &ii, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &jj, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			}
		}
		
	
		/* New clustering */
		Lii=(D[ii][jj]+(S[ii]-S[jj])/(n1-2))/2-LP[ii];
		Ljj=(D[ii][jj]+(S[jj]-S[ii])/(n1-2))/2-LP[jj];

		/* Updating of D */
		if (Lii<0.00001) 
			Lii=0.00005;
		if (Ljj<0.00001)
			Ljj=0.00005;
		L=L+Lii+Ljj;
		LP[ii]=0.5*D[ii][jj];
		Som=Som-(S[ii]+S[jj])/2;
		for (i=1;i<=n1;i++){
			if ((i!=ii)&&(i!=jj)){
				S[i]=S[i]-0.5*(D[i][ii]+D[i][jj]);
				D[i][ii]=(D[i][ii]+D[i][jj])/2;
				D[ii][i]=D[i][ii];
			}
		}
		D[ii][ii]=0;
		S[ii]=0.5*(S[ii]+S[jj])-D[ii][jj];
		if (jj!=n1){
			for (i=1;i<=n1-1;i++){
				D[i][jj]=D[i][n1];
				D[jj][i]=D[n1][i];
			}
			D[jj][jj]=0;
			S[jj]=S[n1];
			LP[jj]=LP[n1];
		}
		/* Updating of DA */
		for (i=1;i<=n;i++){ 
			if (T[i]==ii) 
				T1[i]=T1[i]+Lii;
			if (T[i]==jj)
				T1[i]=T1[i]+Ljj;
		}
		for (j=1;j<=n;j++){ 
			if (T[j]==jj){
				for (i=1;i<=n;i++){
					if (T[i]==ii){
						DA[i][j]=T1[i]+T1[j];
						DA[j][i]=DA[i][j];
					}
				}
			}
		}
		for (j=1;j<=n;j++)
			if (T[j]==jj)  T[j]=ii;

		if (jj!=n1){
			for (j=1;j<=n;j++)
				if (T[j]==n1) T[j]=jj;
		}
	}

	/*Joining the 3 objects remaining */
	l1=(D[1][2]+D[1][3]-D[2][3])/2-LP[1];
	l2=(D[1][2]+D[2][3]-D[1][3])/2-LP[2];
	l3=(D[1][3]+D[2][3]-D[1][2])/2-LP[3];
	if (l1<0.00001) 
		l1=0.00005;
	if (l2<0.00001) 
		l2=0.00005;
	if (l3<0.00001) 
		l3=0.00005;
	L=L+l1+l2+l3;

	for (j=1;j<=n;j++){
		for (i=1;i<=n;i++){
			if ((T[j]==1)&&(T[i]==2)){
				DA[i][j]=T1[i]+T1[j]+l1+l2;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==1)&&(T[i]==3)){
				DA[i][j]=T1[i]+T1[j]+l1+l3;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==2)&&(T[i]==3)){
				DA[i][j]=T1[i]+T1[j]+l2+l3;
				DA[j][i]=DA[i][j];
			}
		}
		DA[j][j]=0;
	}
	
	free(T);
	free(T1);
	free(S);
	free(LP);
	for (i=0;i<=n;i++){ 
  		free(D[i]);   
	} 
	free(D);  
	/*if (myid == 0) {
		endwtime = MPI_Wtime();
		printf("wall clock time = %f\n", endwtime-startwtime);
		totaltime = totaltime + endwtime - startwtime;
	}*/
} /* END NJ */



void construction2(double **DI,double **DA,int *X)
{
	int i,j,k,p,p1,P;
	double *A, *L,B,C,D,a,c,a1,c1,S,S1,S2,S3,G1,M,MR,M1,DIS,DIS1,Su,Sv,Dummy1;
	
	A=(double *) malloc((n+1)*sizeof(double));
	L=(double *) malloc((n+1)*sizeof(double));
	
	Dummy1=100;
	if (DI[X[1]][X[2]]<(n-1)*0.0002){ Dummy1=DI[X[1]][X[2]]; DI[X[1]][X[2]]=(n-1)*0.0002;}
	
	MR=0;
	L[1]=DI[X[1]][X[2]];
	DA[X[1]][X[2]]=L[1];
	DA[X[2]][X[1]]=L[1];
	DA[X[1]][X[1]]=0;
	DA[X[2]][X[2]]=0;
	P=1;
	printf("\navant le for");
	for(k=2;k<=n-1;k++)
	{
		p=1;
		while ((p!=k-1)&&(fabs(DA[X[1]][X[p+1]]+DA[X[1]][X[k]]-DA[X[p+1]][X[k]])<=0.00001))
			p=p+1;
		Su=0;
		Sv=L[1];
		
		for (j=1;j<=P;j++)
		{
			if (j==1)
			{
				B=4*p-2*k;
				C=0;
				D=0;
				G1=0;
				for (i=1;i<=p;i++)
				{
					A[i]=DI[X[i]][X[k+1]]-DA[X[i]][X[k]]+DA[X[1]][X[k]];
					C=C-2*A[i];
					D=D-2*A[i];
					G1=G1+A[i]*A[i];
				}
				
				for (i=p+1;i<=k;i++)
				{
					A[i]=DI[X[i]][X[k+1]]-DA[X[1]][X[i]];
					C=C+2*A[i];
					D=D-2*A[i];
					G1=G1+A[i]*A[i];
				}
			}
			else
			{
				p1=p;
				Su=Sv;
				Sv=Sv+L[j];
				p=p+1;
				
				while (((DA[X[1]][X[p+1]]+DA[X[1]][X[k]]-DA[X[p+1]][X[k]])/2-Sv)<=-0.00001)
				{
				p=p+1;}
				
				S=0;S2=0;S3=0;
				for (i=p1+1;i<=p;i++)
				{
					S=S+A[i];
					S1=DI[X[i]][X[k+1]]-DA[X[i]][X[k]]+DA[X[1]][X[k]];
					S2=S2-A[i]*A[i]+S1*S1;
					A[i]=S1;
					S3=S3+A[i];
				}
				B=4*p-2*k;
				C=C+2*(-S-S3);
				D=D+2*(S-S3);
				G1=G1+S2;
			}
			
			a=Sv;c=0;M=k*a*a+C*a;
			
			if (((c1=-(B*Sv+D)/(2*k))>=-0.00001)&&((M1=k*(Sv*Sv+c1*c1)+B*Sv*c1+C*Sv+D*c1)<M))
			{
				a=Sv;c=c1,M=M1;
			}
			if (((a1=-C/(2*k))>=Su-0.00001)&&(a1<=Sv+0.00001)&&((M1=k*a1*a1+C*a1)<M))
			{
				a=a1;c=0,M=M1;
			}
			if ((2*k!=B)&&(((a1=(B*D-2*k*C)/(4*k*k-B*B))>=Su-0.00001)&&(a1<=Sv+0.00001)&&((c1=(B*C-2*k*D)/(4*k*k-B*B))>=-0.00001)&&((M1=k*(a1*a1+c1*c1)+B*a1*c1+C*a1+D*c1)<M)))
			{
				a=a1;c=c1,M=M1;
			}
			if ((M1=k*Su*Su+C*Su)<M)
			{
				a=Su;c=0,M=M1;
			}
			if (((c1=-(B*Sv+D)/(2*k))>=-0.00001)&&((M1=k*(Su*Su+c1*c1)+B*Su*c1+C*Su+D*c1)<M))
			{
				a=Su;c=c1,M=M1;
			}
			
			
			if ((j==1)||(MR>M+G1))
			{
				MR=M+G1;
				DIS=a;
				DIS1=c;
				
				if (c<=(n-k+1.0)*0.00002) DIS1=(n-k+1.0)*0.00002;
				if (fabs(Sv-Su)<=2*(n-k+1)*0.00002) DIS=Sv-(Sv-Su)/(n-k+1.0);
				else 
				{   
					if (fabs(a-Su)<=(n-k+1.0)*0.00002) DIS=Su+(n-k+1.0)*0.00002;
					if (fabs(a-Sv)<=(n-k+1.0)*0.00002) DIS=Sv-(n-k+1.0)*0.00002;
				} 
			}
		}
		
		i=0;
		S=0;
		while (S<DIS-0.00001)
		{
			i=i+1;
			S=S+L[i];
		}
		
		if ((S==DIS)&&(DIS1>0))
		{
			L[i+1]=DIS1;
			P=i+1;
		}
		
		if ((S>DIS)&&(DIS1>0))
		{
			L[i]=L[i]+DIS-S;
			L[i+1]=DIS1;
			P=i+1;
		}
		
		DA[X[1]][X[k+1]]=DIS+DIS1;
		DA[X[k+1]][X[1]]=DIS+DIS1;
		
		
		
		/*Induction formula */
		
		for (j=2;j<=k;j++)
		{
			if (((DA[X[1]][X[j]]+DA[X[1]][X[k]]-DA[X[k]][X[j]])/2)<=DIS+0.00001)
				DA[X[j]][X[k+1]]=DIS+DIS1+DA[X[j]][X[k]]-DA[X[1]][X[k]];
			
			else
				DA[X[j]][X[k+1]]=DIS1-DIS+DA[X[1]][X[j]];
			DA[X[k+1]][X[j]]=DA[X[j]][X[k+1]];
			DA[X[k+1]][X[k+1]]=0;
		}
		
	}
	printf("\napres le for");
	if (Dummy1<(n-1)*0.0002) DI[X[1]][X[2]]=Dummy1;
	
	free(A);
	free(L);
} /* construction2 */


/* Compute Robinson and Foulds topological distance between two trees 
   given by their distance matrices D and D1*/

int comparb (double **D,double **D1)
{
	struct MATR { unsigned int V:1; };
	struct MATR **A,**A1;
	int i,j,k,p,p1,*P,*P1,*X,*X1,R;
	double S,DIS,DIS1,*L,*L1;
	
	P=(int *)malloc((2*n-2)*sizeof(int));
	P1=(int *)malloc((2*n-2)*sizeof(int));
	X=(int *)malloc((2*n-2)*sizeof(int));
	X1=(int *)malloc((2*n-2)*sizeof(int));
	L=(double *)malloc((2*n-2)*sizeof(double));
	L1=(double *)malloc((2*n-2)*sizeof(double));
        
	A=(struct MATR **)malloc((2*n-2)*sizeof(struct MATR*));
	A1=(struct MATR **)malloc((2*n-2)*sizeof(struct MATR*));      
	
	for (i=0;i<=2*n-3;i++)
	{
		A[i]=(struct MATR *)malloc((n+1)*sizeof(struct MATR));
		A1[i]=(struct MATR *)malloc((n+1)*sizeof(struct MATR));
		
		
		if ((A[i]==NULL)||(A1[i]==NULL)) 
		{printf("probleme de memoire");
			exit(5);
		}
	}
	
	i=1; j=n;
	odp(D,X,&i,&j);
	odp(D1,X1,&i,&j);
	
	for (i=1;i<=2*n-3;i++)
	{
		for (j=0;j<=n;j++)
		{
			A[i][j].V=0;
			A1[i][j].V=0;
		}
	}
	A[1][X[2]].V=1;
	A1[1][X1[2]].V=1;
	P[1]=1;
	P1[1]=1;
	L[1]=D[X[1]][X[2]];
	L1[1]=D1[X1[1]][X1[2]];
	p=1;
	p1=1;
	
	for(k=2;k<=n-1;k++)
	{
		S=0;
		i=0;
		DIS=(D[X[1]][X[k]]+D[X[1]][X[k+1]]-D[X[k]][X[k+1]])/2;
		DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;
		while (S<DIS-0.001)
		{
			i=i+1;
			S=S+L[P[i]];
			A[P[i]][X[k+1]].V=1;
		}
		if ((S-DIS)>=0.001)
		{
			p=p+1;
			for (j=1;j<=k;j++)
				A[p][X[j]].V=A[P[i]][X[j]].V;
			L[P[i]]=L[P[i]]+DIS-S;
			L[p]=S-DIS;
		}
		if (DIS1>=0.001)
		{
			p=p+1;
			A[p][X[k+1]].V=1;
			P[i+1]=p;
			L[p]=DIS1;
		}
		
		S=0;
		i=0;
		DIS=(D1[X1[1]][X1[k]]+D1[X1[1]][X1[k+1]]-D1[X1[k]][X1[k+1]])/2;
		DIS1=(D1[X1[1]][X1[k+1]]+D1[X1[k]][X1[k+1]]-D1[X1[1]][X1[k]])/2;
		while (S<DIS-0.001)
		{
			i=i+1;
			S=S+L1[P1[i]];
			A1[P1[i]][X1[k+1]].V=1;
		}
		
		if ((S-DIS)>=0.001)
		{
			p1=p1+1;
			for (j=1;j<=k;j++)
				A1[p1][X1[j]].V=A1[P1[i]][X1[j]].V;
			L1[P1[i]]=L1[P1[i]]+DIS-S;
			L1[p1]=S-DIS;
		}
		if (DIS1>=0.001)
		{
			p1=p1+1;
			A1[p1][X1[k+1]].V=1;
			P1[i+1]=p1;
			L1[p1]=DIS1;
		}
	}
	
	
	for (i=1;i<=p;i++)
	{
		for (k=1;k<=p1;k++)
		{
			if (A1[k][0].V==0)
			{
				S=0;
				for (j=1;j<=n;j++)
				{
					if (A[i][j].V==A1[k][j].V)
						S=S+1;
				}
				if ((S==0)||(S==n))
				{
					A1[k][0].V=1;
					A[i][0].V=1;
					k=p1;
				}
			}
		}
	}
	R=0;
	for (i=1;i<=p;i++)
	{
		if (A[i][0].V==0)
			R=R+1;
	}
	
	for (i=1;i<=p1;i++)
	{
		if (A1[i][0].V==0)
			R=R+1;
	}
	
	
	free(P);
	free(X);
	free(L);
	free(P1);
	free(X1);
	free(L1);
	
	
	for (i=0;i<=2*n-3;i++)
	{ 
		free(A[i]);
		free(A1[i]);  
	}     
	free(A);
	free(A1);
	
	return R;
} /* END comparb */

/* Compare if two trees given by their distance matrices D and D1 are identical*/

int comparb1 (double **D,double **D1)
{
	int i,j,k,p,p1,*X,*Y;
	double S,S1,DIS,DIS1,DIS0,DIS10,*L,*L1;
	
	
	X=(int *)malloc((n+1)*sizeof(int));
	Y=(int *)malloc((n+1)*sizeof(int));
	L=(double *)malloc((n+1)*sizeof(double));
	L1=(double *)malloc((n+1)*sizeof(double));
	
	
	i=1; j=n;
	odp(D,X,&i,&j);
	
	for (j=1;j<=n;j++)
		Y[j]=0;
	Y[X[1]]=1;
	Y[X[n]]=1;
	for (i=0;i<=n-3;i++)
	{
		S=D1[X[n-i]][X[n-i-1]]-D1[X[n-i-1]][X[1]];
		Y[X[n-i-1]]=1;
		for (j=2;j<=n-i-2;j++)
		{
			S1=D1[X[n-i]][X[j]]-D1[X[j]][X[1]];
			if (S1<S-0.01) 
			{      
				free(Y);
				free(L);
				free(L1); 
				free(X);
				return 0;
			}
		}
	}
	
	L[1]=D[X[1]][X[2]];
	L1[1]=D1[X[1]][X[2]];
	p=1;p1=1;
	for(k=2;k<=n-1;k++)
	{
		
		DIS=(D[X[1]][X[k]]+D[X[1]][X[k+1]]-D[X[k]][X[k+1]])/2;
		DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;
		
		S=0;
		i=0;
		while (S<DIS-0.001)
		{
			i=i+1;
			S=S+L[i];
		}
		
		DIS0=(D1[X[1]][X[k]]+D1[X[1]][X[k+1]]-D1[X[k]][X[k+1]])/2;
		DIS10=(D1[X[1]][X[k+1]]+D1[X[k]][X[k+1]]-D1[X[1]][X[k]])/2;
		
		S1=0;
		j=0;
		while (S1<DIS0-0.001)
		{
			j=j+1;
			S1=S1+L1[j];
		}
		if ((i!=j)||((S-DIS<0.001)&&(S1-DIS0>0.001))||((S-DIS>0.001)&&(S1-DIS0<0.001)))
			
		{  
			free(Y);
			free(L);
			free(L1); 
			free(X);
			return 0;
		}
		
		p=i;
		if ((S-DIS)>=0.001)
		{
			L[i]=L[i]+DIS-S;
		}
		if (DIS1>=0.001)
		{
			p=p+1;
			L[p]=DIS1;
		}
		
		p1=j;
		if ((S1-DIS0)>=0.001)
		{
			L1[j]=L1[j]+DIS0-S1;
		}
		if (DIS10>=0.001)
		{
			p1=p1+1;
			L1[p1]=DIS10;
		}
	}   
	free(Y);
	free(L);
	free(L1); 
	free(X);
	
	return 1;
} /* END comparb1 */

/* Coding the tree metric DA (n,n) by the vectors C (2n-3) and X(n)*/

void coder(double **DA,int *X,double *C)
{
	int i,j;
	
	i=1; j=n;
	odp(DA,X,&i,&j);  
	for (i=1;i<=n-2; i++)
	{
		C[2*i-1]=DA[X[i]][X[i+1]];
		C[2*i]=DA[X[1]][X[i+2]];
		
	}
	C[2*n-3]=DA[X[n-1]][X[n]];
	
} /* END coder */

/* Decoding the tree metric DA (n,n) from the vectors C (2n-3) and X(n)*/

void decoder(double **DA,int *X,double *C)
{
	int i,j;
	
	for (i=1;i<=n-2; i++)
	{
		DA[X[i]][X[i+1]]=C[2*i-1];
		DA[X[1]][X[i+2]]=C[2*i];
		
		DA[X[i+1]][X[i]]=C[2*i-1];
		DA[X[i+2]][X[1]]=C[2*i];
	}
	DA[X[n-1]][X[n]]=C[2*n-3];
	DA[X[n]][X[n-1]]=C[2*n-3];
	
	DA[X[1]][X[1]]=0;
	DA[X[2]][X[2]]=0;
	DA[X[3]][X[3]]=0;
	
	for (i=4;i<=n;i++)
	{
		DA[X[i]][X[i]]=0;
		for (j=2;j<=i-2;j++)
		{
			if (DA[X[1]][X[j]]+DA[X[i-1]][X[i]]<DA[X[1]][X[i]]+DA[X[j]][X[i-1]])
				DA[X[j]][X[i]]=DA[X[1]][X[i]]+DA[X[j]][X[i-1]]-DA[X[1]][X[i-1]];
			else
				DA[X[j]][X[i]]=DA[X[1]][X[j]]+DA[X[i-1]][X[i]]-DA[X[1]][X[i-1]];
			DA[X[i]][X[j]]=DA[X[j]][X[i]];
		}
	}
} /* END decoder */

/* Compute tree statistics and save them in the Output file */
void compute_criteres(double **D,double **DA,double *R,double **W, int optionTR, char **Names)
{
	int i1,j1,*X;
	double s;
	float c;
	char noms[50]; int p;
	char debut[] = "et = ";
	char virgule[] = ",";
	int i,j;
	
	if (optionTR!=3) {
		
		X=(int *)malloc((n+1)*sizeof(int));
		R[1]=0;
		R[2]=0;
		R[3]=W[1][2]*fabs(D[1][2]-DA[1][2]);
		for (i1=1;i1<=n-1;i1++)
		{
			for (j1=i1+1;j1<=n;j1++)
			{
				
				s=fabs(D[i1][j1]-DA[i1][j1]);
				
				R[1]=R[1]+W[i1][j1]*s*s;
				
				if ((W[0][0]<1)||(W[0][0]==2))
				{
					R[2]=R[2]+W[i1][j1]*s;
					if (W[i1][j1]*s>R[3]) R[3]=W[i1][j1]*s;
				}
				
			}
		} 
		if ((W[0][0]<1)||(W[0][0]==2))
		{
			R[2]=2*R[2]/n/(n-1);
			i1=1; j1=n;
			odp(DA,X,&i1,&j1);
			R[0]=DA[X[1]][X[n]];;
			for (i1=1;i1<=n-1;i1++)
				R[0]=R[0]+DA[X[i1]][X[i1+1]];
		}
		R[0]=R[0]/2;
	}
	
	/* Print tree metric matrix if needed*/
	if (W[0][0]==0.1)
	{
		if (optionTR!=2) fprintf(Output1,"\n\nTREE METRIC (ADDITIVE DISTANCE) MATRIX (AD)\n\n");
		else fprintf(Output1,"\n\nRETICULOGRAM DISTANCE MATRIX (AD)\n\n");
		
		
		fprintf(TreeFile,"et = ");
		
		for ( i=1; i<=n; i++)
		{ 	
			
			fprintf(Output1,"%s ",Names[i-1]);
			
			p=0;
			while((Names[i-1][p] !=' ')&&(p<=30)){
				noms[p] = Names[i-1][p];
				p++;
			}
			
			noms[p] = '\0';
			
			if(i!=1) fprintf(TreeFile,",");
			fprintf(TreeFile,"%s",noms);
			
			for ( j=1; j<=n; j++)
			{
				c=(float) DA[i][j];
				fprintf(Output1,"%f ",c);  		
			}
			fprintf(Output1,"\n");    
		}  
	}
	
	if (W[0][0]<1) 
		if (optionTR!=2) fprintf(Output1,"\n\nTHE FOLLOWING STATISTICS ARE AVAILABLE FOR \nA GIVEN DISSIMILARITY (D) AND AN OBTAINED TREE METRIC (AD)\n\n");
	else fprintf(Output1,"\n\nTHE FOLLOWING STATISTICS ARE AVAILABLE FOR \nA GIVEN DISSIMILARITY (D) AND AN OBTAINED RETICULOGRAM DISTANCE (AD)\n\n");
	
	if (optionTR==3) fprintf(Output1,"Only existing values of D were taken into account when computing these statistics\n\n"); 
	
	if ((W[0][0]<1)&&(method==5))
	{
		fprintf(Output1,"\n");
		if (fabs(R[1])<1) fprintf(Output1, "Weighted least-squares coefficient Sum Wij(Dij-ADij)^2  = %12.10f\n",R[1]);   														 
		else fprintf(Output1, "Weighted least-squares coefficient Sum Wij(Dij-ADij)^2  = %f\n",R[1]);   														 
   		fprintf(Output1,"                                    i<j \n");
		if (fabs(R[2])<1) fprintf(Output1,"Average absolute weighted difference  Sum |Wij(Dij-ADij)|/(n(n-1)/2)  = %12.10f\n",R[2]);
		else fprintf(Output1,"Average absolute weighted difference  Sum |Wij(Dij-ADij)|/(n(n-1)/2)  = %f\n",R[2]);
   		fprintf(Output1,"                                      i<j \n");
		if (fabs(R[3])<1) fprintf(Output1,"Maximum weighted absolut difference   Max|Wij(Dij-ADij)| = %12.10f\n",R[3]);
		else			     fprintf(Output1,"Maximum weighted absolut difference   Max|Wij(Dij-ADij)| = %f\n",R[3]);
		fprintf(Output1,"                                      i,j \n");
		if (optionTR!=2) 
		{
			if (fabs(R[0])<1) fprintf(Output1,"Total length of the tree   L  = %12.10f\n",R[0]);
			else              fprintf(Output1,"Total length of the tree   L  = %f\n",R[0]);
		}
	}
	else if ((W[0][0]<1)&&(method<5))
	{
		fprintf(Output1,"\n");
		if (fabs(R[1])<1) fprintf(Output1, "Least-squares coefficient  Sum (Dij-ADij)^2  = %12.10f\n",R[1]);   														 
		else 			 fprintf(Output1, "Least-squares coefficient  Sum (Dij-ADij)^2  = %f\n",R[1]);   														 
		fprintf(Output1,"                            i<j \n");
		if (fabs(R[2])<1) fprintf(Output1,"Average absolute difference  Sum |Dij-ADij|/(n(n-1)/2)  = %12.10f\n",R[2]);
		else 			 fprintf(Output1,"Average absolute difference  Sum |Dij-ADij|/(n(n-1)/2)  = %f\n",R[2]);
		fprintf(Output1,"                             i<j \n");
		if (fabs(R[3])<1) fprintf(Output1,"Maximum absolute difference   Max|Dij-ADij| = %12.10f\n",R[3]);
		else				 fprintf(Output1,"Maximum absolute difference   Max|Dij-ADij| = %f\n",R[3]);
		fprintf(Output1,"                              i,j \n");
		if (optionTR!=2) 
		{
			if (fabs(R[0])<1) fprintf(Output1,"Total length of the tree   L  = %12.10f\n",R[0]);
			else              fprintf(Output1,"Total length of the tree   L  = %f\n",R[0]);
		}
	}
	
	if (optionTR!=3) free(X);
} /* END compute_criteres */ 


/* MW general procedure*/

void parcour2(double **DISS,double **W,double **TM,int *Iternum,long int *ARETE, double *LONGUEUR, int unproc)
{
	int i,j,i1,j1,iopt,jopt,*X,Iternumber,Option;
	double **TM1;
	double EQ,EQmin;
	int *Tiopt, *Tjopt; /* Tableau des iopt et jopt necessaire à la parallelisation */
	double *Teqmin;	/* Tableau des Teqmin necessaire à la parallelisation */

	X=(int *) malloc((n+1)*sizeof(int)); 
	TM1=(double **) malloc((n+1)*sizeof(double*));
	/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1 || unproc==1){
		Tiopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Tjopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Teqmin=(double*)malloc((numprocs+1)*sizeof(double));   
	}
	EQmin=(double)LONG_MAX;


	for (i=0;i<=n;i++){
		TM1[i]=(double*)malloc((n+1)*sizeof(double));    
		if (TM1[i]==NULL){
			printf("probleme de memoire");
			exit(5);
		}  
	}    

	Iternumber=*Iternum;
	Option=floor1(DISS[0][0]);
	X[0]=1;
	
	/* Cas d'un processeur */
	if(numprocs==1 || unproc==1){
		for (i=1;i<=n;i++){ 
			for (j=i+1;j<=n;j++)   {
				X[1]=i;
				X[2]=j;          
				construction(DISS,TM1,X,W);     
				approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
				EQ=0.0;
				for (i1=1;i1<=n;i1++){
					for (j1=i1+1;j1<=n;j1++)
					EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
				} 
				if (((i==1)&&(j==2))||(EQ<EQmin)){
					EQmin=EQ;iopt=i;jopt=j;     
				}         
				if (Option==2) break;  
			}
		}
	}
	/* Cas de plusieurs processeurs */
	else{
		for (i=1;i<=n;i++){ 
			if(i%numprocs == myid){
				for (j=i+1;j<=n;j++){
					X[1]=i;
					X[2]=j;          
					construction(DISS,TM1,X,W);     
					approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
					EQ=0.0;
					for (i1=1;i1<=n;i1++){
						for (j1=i1+1;j1<=n;j1++)
						EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
					} 
					if (((i==1)&&(j==2))||(EQ<EQmin)){
						EQmin=EQ;iopt=i;jopt=j;Teqmin[myid]=EQ;Tiopt[myid]=i;Tjopt[myid]=j;            
					}         
					if (Option==2) break;  
				}
			}	
		}
		/* Synchronisation des résultats */
		if(myid>0){
			MPI_Send( &EQmin, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
			MPI_Send( &iopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
			MPI_Send( &jopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
			MPI_Recv( &iopt, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
			MPI_Recv( &jopt, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
			if(myid<numprocs-1){
				MPI_Send( &iopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &jopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			}
		}
		else{
			for(i=1; i<numprocs;i++){
				MPI_Recv( &Teqmin[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tiopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tjopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
			}
			for(i=1; i<numprocs;i++){
				if (((Tiopt[i]==1)&&(Tjopt[i]==2))||(Teqmin[i]<EQmin)){
						EQmin=Teqmin[i];iopt=Tiopt[i];jopt=Tjopt[i];    
					}
			}
			MPI_Send( &iopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			MPI_Send( &jopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
		}
	}
	if (n>10) Iternumber=2*n-3;
	else Iternumber=15;
	X[1]=iopt;
	X[2]=jopt;  
	construction(DISS,TM1,X,W);
	approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);

	free(X);  
	for (i=0;i<=n;i++)   
	free(TM1[i]);   
	free(TM1);
} /* END parcour2 */


/* ADDTREE method  */
void ADDTREE(double **D1,double **DA, int unproc)
{
	double Som,Smin,Sij,L,Lii,Ljj,l1,l2,l3;
	int *T,i,j,ii,jj,n1,i1,j1;
	double **D,*LP,*S,*T1;
	double *Tsmin; /* Tableau des Smin necessaire à la parallelisation */
	int *Tii, *Tjj; /* Tableau des ii et jj necessaire à la parallelisation */


	D=(double **) malloc((n+1)*sizeof(double*));
	T1=(double *) malloc((n+1)*sizeof(double));
	S=(double *) malloc((n+1)*sizeof(double));
	LP=(double *) malloc((n+1)*sizeof(double));
	T=(int *) malloc((n+1)*sizeof(int));
	/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1 || unproc==1){
		Tii=(int *) malloc((numprocs+1)*sizeof(int));
		Tjj=(int *) malloc((numprocs+1)*sizeof(int));
		Tsmin=(double *) malloc((numprocs+1)*sizeof(double));
	}
	Smin=(double)LONG_MAX;


	for (i=0;i<=n;i++)
	{
		D[i]=(double*)malloc((n+1)*sizeof(double));
		if (D[i]==NULL){
			printf("probleme de memoire");
			exit(5);
		}
	}
	L=0;Som=0;
	for (i=1;i<=n;i++){
		S[i]=0; LP[i]=0;
		for (j=1;j<=n;j++){
			D[i][j]=D1[i][j];
			S[i]=S[i]+D[i][j];
		}
		Som=Som+S[i]/2.0;
		T[i]=i;
		T1[i]=0;
	}

	/* Main procedure */

	for (n1=n;n1>3;n1--){
	/* Research of the best pair of objects (i,j) for clustering*/
		Smin=0.0;   
		
		/* Cas d'un processeur */
		if(numprocs==1 || unproc==1){
			for (i=1;i<=n1-1;i++){
				for (j=i+1;j<=n1;j++){
					Sij=0.0;
					for (i1=1;i1<=n1-1;i1++){
						if ((i1!=i)&&(i1!=j)){
							for (j1=i1+1;j1<=n1;j1++){
								if ((j1!=i)&&(j1!=j)){
									if ((D[i][i1]+D[j][j1]-D[i][j]-D[i1][j1]>=0)&&(D[i][j1]+D[j][i1]-D[i][j]-D[i1][j1]>=0))
										Sij=Sij+1.0;
								}
							}
						}
					}
					if (Sij>Smin){
						Smin = Sij;ii=i; jj=j;
					}
				}	
			}
		}
		/* Cas de plusieurs processeurs */
		else{
			for (i=1;i<=n1-1;i++){
				if(i%numprocs==myid){
					for (j=i+1;j<=n1;j++){
						Sij=0.0;
						for (i1=1;i1<=n1-1;i1++){
							if ((i1!=i)&&(i1!=j)){
								for (j1=i1+1;j1<=n1;j1++){
									if ((j1!=i)&&(j1!=j)){
										if ((D[i][i1]+D[j][j1]-D[i][j]-D[i1][j1]>=0)&&(D[i][j1]+D[j][i1]-D[i][j]-D[i1][j1]>=0))
											Sij=Sij+1.0;
									}
								}
							}
						}
						if (Sij>Smin){
							Smin = Sij;ii=i; jj=j;Tsmin[myid] = Sij;Tii[myid]=i; Tjj[myid]=j;
						}
					}	
				}
			}
			/* Synchronisation des résultats */
			if(myid>0){
				MPI_Send( &Smin, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
				MPI_Send( &ii, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
				MPI_Send( &jj, 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
				MPI_Recv( &ii, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &jj, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				if(myid<numprocs-1){
					MPI_Send( &ii, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
					MPI_Send( &jj, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				}
			}
			else{
				for(i=1; i<numprocs;i++){
					MPI_Recv( &Tsmin[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Tii[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Tjj[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
				}
				for(i=1; i<numprocs;i++){
					if( Tsmin[i]==Smin && ( ii>Tii[i] || (ii==Tii[i]&& jj>Tjj[i]) )){
						ii = Tii[i];jj = Tjj[i];Smin = Tsmin[i];
					}
					if(  Tsmin[i]> Smin ){
						ii = Tii[i];jj = Tjj[i];Smin = Tsmin[i];
					}	
				}
				MPI_Send( &ii, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &jj, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			}
		}
		/* New clustering */
		Lii=(D[ii][jj]+(S[ii]-S[jj])/(n1-2.0))/2-LP[ii];
		Ljj=(D[ii][jj]+(S[jj]-S[ii])/(n1-2.0))/2-LP[jj];
		if (Lii<0.00001) Lii=0.00005;
		if (Ljj<0.00001) Ljj=0.00005;
		L=L+Lii+Ljj;
		LP[ii]=0.5*D[ii][jj];

		Som=Som-(S[ii]+S[jj])/2.0;
		for (i=1;i<=n1;i++){
			if ((i!=ii)&&(i!=jj)){
				S[i]=S[i]-0.5*(D[i][ii]+D[i][jj]);
				D[i][ii]=(D[i][ii]+D[i][jj])/2;
				D[ii][i]=D[i][ii];
			}
		}
		D[ii][ii]=0;
		S[ii]=0.5*(S[ii]+S[jj])-D[ii][jj];
		if (jj!=n1){
			for (i=1;i<=n1-1;i++){
				D[i][jj]=D[i][n1];
				D[jj][i]=D[n1][i];
			}
			D[jj][jj]=0;
			S[jj]=S[n1];
			LP[jj]=LP[n1];
		}
		/* Mise a jour de DA */
		for (i=1;i<=n;i++){
			if (T[i]==ii) T1[i]=T1[i]+Lii;
			if (T[i]==jj) T1[i]=T1[i]+Ljj;   
		}
		for (j=1;j<=n;j++){
			if (T[j]==jj){
				for (i=1;i<=n;i++){
					if (T[i]==ii){
						DA[i][j]=T1[i]+T1[j];
						DA[j][i]=DA[i][j];
					}
				}
			}
		}

		for (j=1;j<=n;j++)
			if (T[j]==jj)  T[j]=ii;

		if (jj!=n1){
			for (j=1;j<=n;j++)
				if (T[j]==n1) T[j]=jj;
		}
	}

	/*3 objects remaining*/

	l1=(D[1][2]+D[1][3]-D[2][3])/2.0-LP[1];
	l2=(D[1][2]+D[2][3]-D[1][3])/2.0-LP[2];
	l3=(D[1][3]+D[2][3]-D[1][2])/2.0-LP[3];
	if (l1<0.00001) l1=0.00005;
	if (l2<0.00001) l2=0.00005;
	if (l3<0.00001) l3=0.00005;

	L=L+l1+l2+l3;

	for (j=1;j<=n;j++){
		for (i=1;i<=n;i++){
			if ((T[j]==1)&&(T[i]==2)){
				DA[i][j]=T1[i]+T1[j]+l1+l2;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==1)&&(T[i]==3)){
				DA[i][j]=T1[i]+T1[j]+l1+l3;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==2)&&(T[i]==3)){
				DA[i][j]=T1[i]+T1[j]+l2+l3;
				DA[j][i]=DA[i][j];
			}
		}
		DA[j][j]=0;
	}
	free(T);
	free(T1);
	free(S);
	free(LP);

	for (i=0;i<=n;i++){ 
		free(D[i]);   
	} 
	free(D);
} /* END ADDTREE */


/* Circular order reconstruction, main procedure*/

void parcour1(double **DISS,double **W,double **TM,int *Iternum,long int *ARETE, double *LONGUEUR, int unproc)
{
	int i,j,i1,j1,iopt,jopt,*X,Iternumber;
	double **TM1;
	double EQ,EQmin;
	
	double *Teqmin; /* Tableau des EQmin necessaire à la parallelisation */
	int *Tiopt, *Tjopt; /* Tableau des iopt et jopt necessaire à la parallelisation */


	X=(int *) malloc((n+1)*sizeof(int)); 
	TM1=(double **) malloc((n+1)*sizeof(double*));
	/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1 || unproc==1){
		Tiopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Tjopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Teqmin =(double*)malloc((numprocs+1)*sizeof(double));
	}
	EQmin=(double)LONG_MAX;


	for (i=0;i<=n;i++){
		TM1[i]=(double*)malloc((n+1)*sizeof(double));
		if (TM1[i]==NULL){
			printf("probleme de memoire");
			exit(5);
		}  
	}

	Iternumber=*Iternum;
	/* Cas d'un processeur */
	if(numprocs == 1 || unproc==1){
		for (i=1;i<=n;i++){
			for (j=1;j<=n;j++){ 
				if (i!=j){
					odp(DISS,X,&i,&j);
					construction2(DISS,TM1,X);
					approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
					EQ=0.0;
					for (i1=1;i1<=n;i1++){
						for (j1=i1+1;j1<=n;j1++)
							EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
					} 
					if (((i==1)&&(j==2))||(EQ<EQmin)){
						EQmin=EQ;iopt=i;jopt=j;   
					}
				}
			}
		}
	}
	/* Cas de plusieurs processeurs */
	else{
		for (i=1;i<=n;i++){
			if(i%numprocs == myid){
				for (j=1;j<=n;j++){ 
					if (i!=j){
						odp(DISS,X,&i,&j);
						construction2(DISS,TM1,X);
						approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
						EQ=0.0;
						for (i1=1;i1<=n;i1++){
							for (j1=i1+1;j1<=n;j1++)
								EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
						} 
						if (((i==1)&&(j==2))||(EQ<EQmin)){
							EQmin=EQ;iopt=i;jopt=j;Teqmin[myid]=EQ;Tiopt[myid]=i;Tjopt[myid]=j;         
						}
					}
				}
			}
		}
		/* Synchronisation des resultats */
		if(myid>0){
			MPI_Send( &EQmin, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
			MPI_Send( &iopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
			MPI_Send( &jopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
			MPI_Recv( &iopt, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
			MPI_Recv( &jopt, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
			if(myid<numprocs-1){
				MPI_Send( &iopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &jopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			}
		}
		else{
			for(i=1; i<numprocs;i++){
				MPI_Recv( &Teqmin[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tiopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tjopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
			}
			for(i=1; i<numprocs;i++){
				if (((Tiopt[i]==1)&&(Tjopt[i]==2))||(Teqmin[i]<EQmin)){
						EQmin=Teqmin[i];iopt=Tiopt[i];jopt=Tjopt[i];    
					}
			}
			MPI_Send( &iopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			MPI_Send( &jopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
		}
		
	}
	if (n>10) Iternumber=2*n-3;
	else Iternumber=15;
	odp(DISS,X,&iopt,&jopt);
	construction2(DISS,TM1,X); 
	approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);

	//MPI_Barrier(MPI_COMM_WORLD);
	/*if (myid == 0) {
		endwtime = MPI_Wtime();
		totaltime = totaltime + endwtime - startwtime;
	}*/
	free(X);  
	for (i=0;i<=n;i++)   
		free(TM1[i]);   
	free(TM1);
}

/* Compute the list of tree edges of a given tree distance DI */
void Tree_edges (double **DI, long int *ARETE, double *LONGUEUR)
{
	
	struct EDGE { unsigned int U; unsigned int V; double LN;};
	struct EDGE *Path,*Tree;
	int i,j,k,p,P,*X;
	double S,DIS,DIS1,*L,**D;
	
	X=(int *)malloc((n+1)*sizeof(int));  
	L=(double *)malloc((n+1)*sizeof(double));
	Tree=(struct EDGE *)malloc((2*n-2)*sizeof(struct EDGE));
	Path=(struct EDGE *)malloc((n+2)*sizeof(struct EDGE));
	
	
	D=(double **) malloc((n+1)*sizeof(double*));
	
	for (i=0;i<=n;i++)
	{
		D[i]=(double*)malloc((n+1)*sizeof(double)); 
		
		if (D[i]==NULL)
		{printf("probleme de memoire");
			/* MyDebugStr("\pData matrix is too large ");*/
			exit(5);
		}
	}
	
	
	i=1; j=n;
	odp(DI,X,&i,&j);
	
	for (i=1;i<=n;i++)
	{ 
		for (j=1;j<=n;j++) 
			D[i][j]=DI[i][j];
	}  
	
	L[1]=D[X[1]][X[2]];
	Path[1].U=X[1];
	Path[1].V=X[2];
	
	p=0;
	P=1;
	
	for(k=2;k<=n-1;k++)
	{
		
		DIS=(D[X[1]][X[k]]+D[X[1]][X[k+1]]-D[X[k]][X[k+1]])/2;
		DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;
		
		S=0.0;
		i=0;
		if (DIS>0.00001)
		{
			while (S<DIS-0.00001)
			{
				i=i+1;
				S=S+L[i];
			}
		}
		else { DIS=0; i=1; }
		
		Tree[p+1].U=n+k-1;
		Tree[p+1].V=Path[i].V;
		Tree[p+1].LN=S-DIS;
		if (Tree[p+1].LN<0) Tree[p+1].LN=0; 
		
		for (j=i+1;j<=P;j++)
		{
			Tree[p+j-i+1].U=Path[j].U;
			Tree[p+j-i+1].V=Path[j].V;
			Tree[p+j-i+1].LN=L[j];
			if (L[j]<0) L[j]=0; 
		}
		p=p+P-i+1;
		
		
		
		Path[i].V=n+k-1;
		Path[i+1].U=n+k-1;
		Path[i+1].V=X[k+1];
		L[i]=L[i]+DIS-S;
		L[i+1]=DIS1;
		P=i+1;
	}
	
	for (i=1;i<=P;i++) 
	{
		Tree[p+i].U=Path[i].U;
		Tree[p+i].V=Path[i].V;
		Tree[p+i].LN=L[i];
	}	
	for (i=1;i<=2*n-3;i++)
	{
		if (fabs(Tree[i].LN-0.00001)<=0.00001)
			Tree[i].LN=0.0;
		ARETE[2*i-2]=Tree[i].U;
		ARETE[2*i-1]=Tree[i].V;
		LONGUEUR[i-1]=Tree[i].LN;     
	} 
	
	free(X);
	free(Tree);
	free(L);
	free(Path);
	
	
	for (i=0;i<=n;i++)    
		free(D[i]);
	
	free(D);
}

/* UNJ method*/

void UNJ(double **D1,double **DA, int unproc)
{
	double **D,*T1,*S,*LP,Som,Smin,Sij,L,Lii,Ljj,l1,l2,l3,S1,*NUMBER;
	int *T,i,j,ii,jj,n1;
	
	double *Tsmin; /* Tableau des Smin necessaire à la parallelisation */
	int *Tii, *Tjj; /* Tableau des ii et jj necessaire à la parallelisation */


	D=(double **) malloc((n+1)*sizeof(double*));
	T1=(double *) malloc((n+1)*sizeof(double));
	S=(double *) malloc((n+1)*sizeof(double));
	LP=(double *) malloc((n+1)*sizeof(double));
	T=(int *) malloc((n+1)*sizeof(int));
	NUMBER=(double *) malloc((n+1)*sizeof(double));

		/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1 || unproc==1){
		Tjj=(int *) malloc((numprocs+1)*sizeof(int));
		Tii=(int *) malloc((numprocs+1)*sizeof(int));
		Tsmin=(double *) malloc((numprocs+1)*sizeof(double));
	}
	Smin=(double)LONG_MAX;

	for (i=0;i<=n;i++){
		D[i]=(double*)malloc((n+1)*sizeof(double));
		if (D[i]==NULL){
			printf("probleme de memoire");
			exit(5);
		}
	}
	L=0;
	Som=0;
	for (i=1;i<=n;i++){
		S[i]=0; LP[i]=0; NUMBER[i]=1.0; 
		for (j=1;j<=n;j++){
			D[i][j]=D1[i][j];
			S[i]=S[i]+D[i][j];
		}
		T[i]=i;
		T1[i]=0;
	}

	/* Main procedure */

	for (n1=n;n1>3;n1--){
		/* Research of the best pair of objects (i,j) to be agglomerated  */
		Smin=-S[1]-S[2]+D[1][2]*(n1-2);

		/* Cas d'un seul processeur */
		if(numprocs==1 || unproc==1){
			for (i=1;i<=n1-1;i++){
				for (j=i+1;j<=n1;j++){
					Sij=-S[i]-S[j]+D[i][j]*(n1-2);
					if (Sij<=Smin){
						Smin=Sij;ii=i;jj=j;
					}
				}
			}
		}
		/* Cas de plusieurs processeurs */
		else{
			for (i=1;i<=n1-1;i++){
				if(i%numprocs==myid){
					for (j=i+1;j<=n1;j++){
						Sij=-S[i]-S[j]+D[i][j]*(n1-2);
						if (Sij<=Smin){
							Smin = Sij;ii=i; jj=j;Tsmin[myid] = Sij;Tii[myid]=i; Tjj[myid]=j;
						}
					}
				}
			}
			if(myid>0){
				MPI_Send( &Smin, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
				MPI_Send( &ii, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
				MPI_Send( &jj, 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
				MPI_Recv( &ii, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &jj, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				if(myid<numprocs-1){
					MPI_Send( &ii, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
					MPI_Send( &jj, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				}
			}
			else{
				for(i=1; i<numprocs;i++){
					MPI_Recv( &Tsmin[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Tii[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Tjj[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
				}
				for(i=1; i<numprocs;i++){
					if(Tsmin[i]<=Smin){
						ii = Tii[i];jj = Tjj[i];Smin = Tsmin[i];
					}	
					if( Tsmin[i]==Smin && ( ii<Tii[i] || (ii==Tii[i]&& jj<Tjj[i]) )){
						ii = Tii[i];jj = Tjj[i];Smin = Tsmin[i];
					}
				}
				MPI_Send( &ii, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &jj, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			}
		}
		/* New clustering*/

		S1=0;
		for (i=1;i<=n1;i++){
			if ((i!=ii)&&(i!=jj))
				S1=S1+NUMBER[i]*(D[ii][i]-D[jj][i]);
		}    
		Lii=S1/(2.0*(n-NUMBER[ii]-NUMBER[jj]))+0.5*D[ii][jj];
		Ljj=-S1/(2.0*(n-NUMBER[ii]-NUMBER[jj]))+0.5*D[ii][jj]; 

		if (Lii<0.00001) Lii=0.00005;
		if (Ljj<0.00001) Ljj=0.00005;
		L=L+Lii+Ljj;

		S[ii]=0;
		for (i=1;i<=n1;i++){
			if ((i!=ii)&&(i!=jj)){
				S1=D[i][ii];
				D[i][ii]=(NUMBER[ii]*(D[i][ii]-Lii)+NUMBER[jj]*(D[i][jj]-Ljj))/(NUMBER[ii]+NUMBER[jj]);            
				S[i]=S[i]+D[i][ii]-(D[i][jj]+S1);     
				D[ii][i]=D[i][ii];
				S[ii]=S[ii]+D[i][ii];
			}
		}
		S1=S[ii];
		D[ii][ii]=0;
		NUMBER[ii]=NUMBER[ii]+NUMBER[jj];

		if (jj!=n1){
			for (i=1;i<=n1-1;i++){
				D[i][jj]=D[i][n1];
				D[jj][i]=D[n1][i];
			}
			D[jj][jj]=0;
			S[jj]=S[n1];
			NUMBER[jj]=NUMBER[n1];
		}

		for (i=1;i<=n;i++){
			if (T[i]==ii) T1[i]=T1[i]+Lii;
			if (T[i]==jj) T1[i]=T1[i]+Ljj;
		}

		for (j=1;j<=n;j++){ 
			if (T[j]==jj){
				for (i=1;i<=n;i++){
					if (T[i]==ii){
						DA[i][j]=T1[i]+T1[j];
						DA[j][i]=DA[i][j];
					}
				}
			}
		}

		for (j=1;j<=n;j++)
			if (T[j]==jj)  T[j]=ii;
		if (jj!=n1){
			for (j=1;j<=n;j++)
				if (T[j]==n1) T[j]=jj;
		}
	}

	/*3 objects remaining */

	l1=(D[1][2]+D[1][3]-D[2][3])/2;
	l2=(D[1][2]+D[2][3]-D[1][3])/2;
	l3=(D[1][3]+D[2][3]-D[1][2])/2;
	if (l1<0.00001) l1=0.00005;
	if (l2<0.00001) l2=0.00005;
	if (l3<0.00001) l3=0.00005;
	L=L+l1+l2+l3;

	for (j=1;j<=n;j++){
		for (i=1;i<=n;i++){
			if ((T[j]==1)&&(T[i]==2)){
				DA[i][j]=T1[i]+T1[j]+l1+l2;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==1)&&(T[i]==3)){
				DA[i][j]=T1[i]+T1[j]+l1+l3;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==2)&&(T[i]==3)){
				DA[i][j]=T1[i]+T1[j]+l2+l3;
				DA[j][i]=DA[i][j];
			}
		}
		DA[j][j]=0;
	}
	/*if (myid == 0) {
		endwtime = MPI_Wtime();
		printf("wall clock time = %f\n", endwtime-startwtime);	
	}*/
	free(T);
	free(T1);
	free(S);
	free(LP);
	free(NUMBER);
	for (i=0;i<=n;i++){ 
		free(D[i]);   
	} 
	free(D);
}

/* returns the closest integer value of x*/
int floor1(double x)
{  
	int i;
	
	if (ceil(x)-floor(x)==2.0) i=(int)x; 
	else if (fabs(x-floor(x)) > fabs(x-ceil(x))) i=(int)ceil(x);
	else i=(int)floor(x);
	return i;
} 

/* Procedure for the polishing the tree metric TM subject to the dissimilarity DISS and the weight function W,
   Iternum is the number of iteration in the Gauss-Seidel iterative procedure applyed, 
   TMnew is the matrix of the polished tree metric, 
   ARETE is the vector of tree edges, and  LONGUEUR is the vector of their lengths.*/

void approx_arb2(double **DISS,double **TM,double **TMnew,double **W, int *Iternum, long int *ARETE, double *LONGUEUR)
{
	
	long int  *Level, *Score, *EX1, *EX2, *Succ1, *Succ2, *Succ11, *Succ22, *Neighbour1, *Neighbour2;
	int i,j,k,k1,j1,i1,i2,*Flag, **Part, **Vertices, *VertexNumber;
	
	double *L, **B, *C, Sum, Sum1, *Path, EQ, l;
	
	/*Variable declaration*/
	
	Level=(long int *)malloc((2*n-1)*sizeof(long int));  
	Score=(long int *)malloc((2*n-1)*sizeof(long int));  
	EX1=(long int *)malloc((2*n-1)*sizeof(long int));  
	EX2=(long int *)malloc((2*n-1)*sizeof(long int));
	Succ1=(long int *)malloc((2*n-1)*sizeof(long int)); 
	Succ2=(long int *)malloc((2*n-1)*sizeof(long int));
	Succ11=(long int *)malloc((2*n-1)*sizeof(long int)); 
	Succ22=(long int *)malloc((2*n-1)*sizeof(long int));
	C=(double *)malloc((2*n-1)*sizeof(double));
	B=(double **)malloc((2*n-1)*sizeof(double*));
	
	VertexNumber=(int *)malloc((2*n-1)*sizeof(int));    
	Vertices=(int **)malloc((2*n-1)*sizeof(int*));
	Part=(int **)malloc((2*n-1)*sizeof(int*));
	
	Neighbour1=(long int *)malloc((2*n-1)*sizeof(long int));
	Neighbour2=(long int *)malloc((2*n-1)*sizeof(long int));
	Flag=(int *)malloc((2*n-1)*sizeof(int));
	
	L=(double *)malloc((2*n-1)*sizeof(double));
	Path=(double *)malloc((n+1)*sizeof(double));  
	
	for (i=0;i<=2*n-2;i++)
	{
		Vertices[i]=(int*)malloc((n+1)*sizeof(int));  
		B[i]=(double*)malloc((2*n-1)*sizeof(double));
		Part[i]=(int*)malloc((2*n-1)*sizeof(int));   
		if ((B[i]==NULL)||(Part[i]==NULL)||(Vertices[i]==NULL))
		{
			printf ("Data matrix is too large ");
			exit(5);
		}
	}   
	
	/* Variable initialisation*/
	
	Tree_edges (TM,ARETE,LONGUEUR);
	kt = 0;
        
	for (i=1;i<=2*n-3-kt;i++)
	{        
		EX1[i]=ARETE[2*i-2];
		EX2[i]=ARETE[2*i-1];
		
		L[i]=LONGUEUR[i-1];
		Flag[i]=0;
		
		if (EX1[i]<=n)
			{ Score[EX1[i]]=3; Succ1[EX1[i]]=0; Succ2[EX1[i]]=0; Score[EX2[i]]=0; }
		
		else if (EX2[i]<=n)
			{ Score[EX2[i]]=3; Succ1[EX2[i]]=0; Succ2[EX2[i]]=0; Score[EX1[i]]=0; }
		
		else { Score[EX1[i]]=0; Score[EX2[i]]=0; } 
	}
	for (i=1;i<=2*n-3-kt;i++)
	{ 
		if (i<=n) { Path[i]=0.0; TMnew[i][i]=0; }
		if ((EX1[i]<=n)||(EX2[i]<=n))
		{
			VertexNumber[i]=1;
			if (EX1[i]<=n) Vertices[i][1]=EX1[i];
			else Vertices[i][1]=EX2[i];     
		}
		else VertexNumber[i]=0;
		
		for (j=1;j<=2*n-3-kt;j++)
		{
			Part[i][j]=0;
			Part[j][i]=0;   
		} 
	}
	
	
	/*   Filling in the vectors Level(2n-3), VertexNumber(2n-3),
	     and the matrices Part(2n-3, 2n-3), Vertices(2n-3,n)   */ 
	
	j=1; 
	while (j<2*n-2-kt)
	{    
		for (i=1; i<=2*n-3-kt; i++)
		{
			
			if (Score[EX1[i]]==5) Score[EX1[i]]=3;
			if (Score[EX2[i]]==5) Score[EX2[i]]=3;
			
			if (((Score[EX1[i]]==3)&&(Score[EX2[i]]!=4))||((Score[EX2[i]]==3)&&(Score[EX1[i]]!=4))) 
			{ 
				if (Score[EX1[i]]==3) { k=EX2[i]; k1=EX1[i]; }
				else { k=EX1[i]; k1=EX2[i]; }
				
				if (Score[k]<2) 
				{
					Score[k]=Score[k]+1;
					if (Score[k]==1) { Succ1[k]=i; Neighbour1[k]=k1; Neighbour2[k]=0;}
					else  { Succ2[k]=i; Neighbour2[k]=k1; }
				}
			}
			
		}
		
		for (i=1; i<=2*n-3-kt; i++)
		{
			if (((Score[EX1[i]]==1)&&(Score[EX2[i]]==3))||((Score[EX1[i]]==2)&&(Score[EX2[i]]==3))||
				((Score[EX1[i]]==3)&&(Score[EX2[i]]==1))||((Score[EX1[i]]==3)&&(Score[EX2[i]]==2)))
			
			{ 
				
				if ((Score[EX1[i]]==1)||(Score[EX1[i]]==2)) { k=EX1[i]; k1=EX2[i]; }
				else { k=EX2[i]; k1=EX1[i]; } 
				
				if(Flag[i]==0) 
				{ 
					Succ11[i]=Succ1[k1];
					Succ22[i]=Succ2[k1];
					Level[j]=i; j=j+1; Score[k1]=4;
					Flag[i]=1;
					if ((Score[Neighbour1[k]]==4)&&(Score[Neighbour2[k]]==4)) Score[k]=5;
					
					if (j>n+1)
					{              
						VertexNumber[i]=VertexNumber[Succ11[i]]+VertexNumber[Succ22[i]];
						for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)
							Vertices[i][i1]=Vertices[Succ11[i]][i1];
						
						for (i1=VertexNumber[Succ11[i]]+1; i1<=VertexNumber[i]; i1++)
							Vertices[i][i1]=Vertices[Succ22[i]][i1-VertexNumber[Succ11[i]]];
						
						
						Part[Succ11[i]][i]=1; Part[Succ22[i]][i]=1;
						Part[i][Succ11[i]]=1; Part[i][Succ22[i]]=1; 
						for (i1=1; i1<=2*n-3; i1++)
						{          
							if ((Part[Succ11[i]][i1]==1)||(Part[Succ22[i]][i1]==1))
							{
								Part[i][i1]=1; Part[i1][i]=1;
							}
						}
					} 
					
					
				}
				
				if(Score[k]==2) 
				{ 
					Score[k]=5;
					Succ11[Succ2[k]]=Succ1[Neighbour2[k]];
					Succ22[Succ2[k]]=Succ2[Neighbour2[k]];
					Level[j]=Succ2[k]; j=j+1; Score[Neighbour2[k]]=4;
					Flag[Succ2[k]]=1;       
					
					j1=i;
					i=Succ2[k];
					
					if (j>n+1)
					{              
						VertexNumber[i]=VertexNumber[Succ11[i]]+VertexNumber[Succ22[i]];
						for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)
							Vertices[i][i1]=Vertices[Succ11[i]][i1];
						
						for (i1=VertexNumber[Succ11[i]]+1; i1<=VertexNumber[i]; i1++)
							Vertices[i][i1]=Vertices[Succ22[i]][i1-VertexNumber[Succ11[i]]];
						
						
						Part[Succ11[i]][i]=1; Part[Succ22[i]][i]=1;
						Part[i][Succ11[i]]=1; Part[i][Succ22[i]]=1; 
						for (i1=1; i1<=2*n-3; i1++)
						{          
							if ((Part[Succ11[i]][i1]==1)||(Part[Succ22[i]][i1]==1))
							{
								Part[i][i1]=1; Part[i1][i]=1;
							}
						}
					} 
					i=j1;            
					
				} 
				
			} 
			
			
		} 
		
		for (i=1; i<=2*n-3; i++)
		{ 
			if (((Score[EX1[i]]==5)&&(Score[EX2[i]]==3))||((Score[EX1[i]]==3)&&(Score[EX2[i]]==5))||
				((Score[EX1[i]]==5)&&(Score[EX2[i]]==5)))
			{      
				if ((Score[EX1[i]]==3)||((Score[EX1[i]]==5)&&(Score[EX2[i]]==5)))
					{ k=EX1[i]; k1=EX2[i]; }
				else { k=EX2[i]; k1=EX1[i]; } 
				
				Level[j]=i; j=j+1;    
				
				Succ11[i]=Succ1[k];
				Succ22[i]=Succ2[k];
				
				Score[k]=4;
				Score[k1]=4;
				
				if (j>n+1)
				{              
					VertexNumber[i]=VertexNumber[Succ11[i]]+VertexNumber[Succ22[i]];
					for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)
						Vertices[i][i1]=Vertices[Succ11[i]][i1];
					
					for (i1=VertexNumber[Succ11[i]]+1; i1<=VertexNumber[i]; i1++)
						Vertices[i][i1]=Vertices[Succ22[i]][i1-VertexNumber[Succ11[i]]];
					
					
					Part[Succ11[i]][i]=1; Part[Succ22[i]][i]=1;
					Part[i][Succ11[i]]=1; Part[i][Succ22[i]]=1; 
					for (i1=1; i1<=2*n-3; i1++)
					{          
						if ((Part[Succ11[i]][i1]==1)||(Part[Succ22[i]][i1]==1))
						{
							Part[i][i1]=1; Part[i1][i]=1;
						}
					}
				} 
			}       
		}  
	}
	
	
	
	/* Filling in the matrix B(2n-3,2n-3) and the vector C(2n-3)  */
	
	for (i=1; i<=2*n-3; i++)
	{
		B[i][i]=0;
		C[i]=0;
		for (j=1; j<=2*n-3; j++)
		{   
			if ((EX1[i]<=n)||(EX2[i]<=n))
			{
				if (EX1[i]<=n) k=EX1[i];
				else k=EX2[i];
				if (((EX1[j]<=n)||(EX2[j]<=n))&&(i!=j))       
				{           
					if (EX1[j]<=n) k1=EX1[j];
					else k1=EX2[j];
					
					B[i][j]=W[k][k1];
					B[j][i]=W[k][k1];
					B[i][i]=B[i][i]+B[i][j];
					C[i]=C[i]+DISS[k][k1]*W[k][k1];     
				}
			}
		}
	}
	
	for (k=n+1; k<=2*n-3; k++)
	{   
		i=Level[k];
		Sum=0.0;
		Sum1=0.0;
		for (j=1; j<=VertexNumber[Succ11[i]]; j++) 
		{
			for (j1=1; j1<=VertexNumber[Succ22[i]]; j1++)     
			{
				Sum=Sum+2*W[Vertices[Succ11[i]][j]][Vertices[Succ22[i]][j1]];
				Sum1=Sum1+2*DISS[Vertices[Succ11[i]][j]][Vertices[Succ22[i]][j1]]*
				W[Vertices[Succ11[i]][j]][Vertices[Succ22[i]][j1]];
			}
		}     
		B[i][i]=B[Succ11[i]][Succ11[i]]+B[Succ22[i]][Succ22[i]]-Sum;
		C[i]=C[Succ11[i]]+C[Succ22[i]]-Sum1;  
	}
	
	
	for (j1=n+1; j1<=2*n-3; j1++)
	{   
		i=Level[j1];
		for (i1=1; i1<j1; i1++)
		{       
			j=Level[i1];
			if ((j==Succ11[i])||(j==Succ22[i]))
			{  
				if (j==Succ11[i]) k=Succ22[i];
				else k=Succ11[i]; 
				B[i][j]=(B[i][i]+B[j][j]-B[k][k])/2;
			}  
			
			else if (Part[i][j]==1)
			{        
				if (Part[j][Succ11[i]]==1)
					B[i][j]=B[j][Succ11[i]]-B[j][Succ22[i]];
				else 
					B[i][j]=B[j][Succ22[i]]-B[j][Succ11[i]];    
			}
			
			else    
			{        
				B[i][j]=B[Succ11[i]][j]+B[Succ22[i]][j];        
			}
			B[j][i]=B[i][j];
			
		}
	}  
	
	
	/* Gausse-Seidel iteretive procedure to polish vector of lengths L(2n-3)*/
	
	for(k=1;k<=*Iternum;k++)
	{    
		
		
		for(i=1;i<=2*n-3;i++)
		{
			EQ=0;
			Sum=0;
			Sum1=0;
			for (j=1;j<=2*n-3;j++)
			{
				if (j>=i+1) Sum=Sum+B[i][j]*L[j];
				if (j<=i-1) Sum1=Sum1+B[i][j]*L[j];
			}
			l=L[i];
			L[i]=(-Sum-Sum1+C[i])/B[i][i];
			if (L[i]<0) L[i]=0.0;
			EQ=EQ+sqrt((L[i]-l)*(L[i]-l));   
		}    
	}
	
	for (i=1;i<=2*n-3;i++)
	{
		LONGUEUR[i-1]=L[i];
	}
	
	/* Computation of the new tree distance matrix TMnew(n,n) from the list of edges*/
	
	for (j=n+1; j<=2*n-3; j++)
	{   
		i=Level[j];
		
		for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)   
			Path[Vertices[Succ11[i]][i1]]=Path[Vertices[Succ11[i]][i1]]+L[Succ11[i]];  
		
		for (i1=1; i1<=VertexNumber[Succ22[i]]; i1++)     
			Path[Vertices[Succ22[i]][i1]]=Path[Vertices[Succ22[i]][i1]]+L[Succ22[i]];   
		
		for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)   
		{   
			k=Vertices[Succ11[i]][i1];
			for (j1=1; j1<=VertexNumber[Succ22[i]]; j1++)     
			{            
				k1=Vertices[Succ22[i]][j1];
				TMnew[k][k1]=Path[k]+Path[k1];
				TMnew[k1][k]=TMnew[k][k1];
			}
		}              
	}
	
	i=Level[2*n-3];
	if ((EX1[i]==EX1[Succ11[i]])||(EX1[i]==EX2[Succ11[i]])) k=EX2[i];
	else k=EX1[i];
	
	j1=1;
	for (i1=1; i1<=2*n-4; i1++)
	{    
		k1=Level[i1];
		if ((EX1[k1]==k)||(EX2[k1]==k))
		{
			if (j1==1) { i=k1; j1=2; } 
			else { j=k1; break; }  
		}
	}
	
	for (i1=1; i1<=VertexNumber[i]; i1++)   
	{    
		for (j1=1; j1<=VertexNumber[j]; j1++)     
		{   
			k=Vertices[i][i1];
			k1=Vertices[j][j1];
			TMnew[k][k1]=Path[k]+Path[k1]+L[i]+L[j];
			TMnew[k1][k]=TMnew[k][k1];
		}
	}        
	
	i2=Level[2*n-3];
	for (i1=1; i1<=VertexNumber[i2]; i1++)   
	{             
		k=Vertices[i2][i1];
		for (j1=1; j1<=VertexNumber[i]; j1++)     
		{   
			k1=Vertices[i][j1];
			TMnew[k][k1]=Path[k]+Path[k1]+L[i2]+L[i];
			TMnew[k1][k]=TMnew[k][k1];
		}     
		for (j1=1; j1<=VertexNumber[j]; j1++)     
		{   
			k1=Vertices[j][j1];
			TMnew[k][k1]=Path[k]+Path[k1]+L[i2]+L[j];
			TMnew[k1][k]=TMnew[k][k1];
		}          
	}   
	
	
	free(Level);  
	free(Score);  
	free(EX1);  
	free(EX2);
	free(Succ1); 
	free(Succ2);
	free(Succ11); 
	free(Succ22);
	free(C);  
	free(VertexNumber);     
	free(Neighbour1);
	free(Neighbour2);
	free(Flag);  
	free(L);
	free(Path);
	
	for (i=0;i<=2*n-2;i++)
	{ 
		free(B[i]);
		free(Part[i]);
		free(Vertices[i]);  	   
	} 
	free(B);
	free(Part);
	free(Vertices); 
}
/* Scaling of the dissimilarity entries to avoid too small or too large values*/
void Scaling(double **DI, double *power)
{
	int i,j;
	double MAX, MIN, POWER;
	
	MAX=fabs(DI[1][2]);
	MIN=fabs(DI[1][2]);
	for (i=1;i<=n;i++)
	{
		for (j=i+1;j<=n;j++)
		{
			if (fabs(DI[i][j])>MAX) MAX=fabs(DI[i][j]);
			if (fabs(DI[i][j])<MIN) MIN=fabs(DI[i][j]);
		} 
	}   	
	
	POWER=0.0;
	while (MAX<=pow(10.0,-POWER))
		POWER=POWER+1;    
	
	if (POWER!=0)
	{
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{  
				DI[i][j]=DI[i][j]*pow(10.0,POWER);
				DI[j][i]=DI[i][j];
			}
		}
	}
        
	if ((POWER==0)&&(MAX>100000.0))
	{
		POWER=-2;
		while (MAX>pow(10.0,-POWER))
			POWER=POWER-1;  
		if (POWER!=0)
		{
			for (i=1;i<=n;i++)
			{
				for (j=i+1;j<=n;j++)
				{  
					DI[i][j]=DI[i][j]*pow(10.0,POWER);
					DI[j][i]=DI[i][j];
				}
			}   
		}    
	}
	*power=POWER;
}

/* The procedure inverse to Scaling to restore true dissimilarity and tree metric values */
void Scaling1(double **DI,double **DA,double *RAJ,double *LONGUEUR, double *power)
{
	int i,j;
	double POWER;
	
	if (*power!=0)
	{
		POWER=*power;
		
		RAJ[1]=RAJ[1]*pow(10.0,-2*POWER);
		RAJ[2]=RAJ[2]*pow(10.0,-POWER);
		RAJ[3]=RAJ[3]*pow(10.0,-POWER);
		RAJ[0]=RAJ[0]*pow(10.0,-POWER);
		
		for (i=1;i<=n;i++)
		{         
			for (j=i+1;j<=n;j++)
			{  
				DI[i][j]=DI[i][j]*pow(10.0,-POWER);
				DI[j][i]=DI[i][j];
				DA[i][j]=DA[i][j]*pow(10.0,-POWER);
				DA[j][i]=DA[i][j];
			}
		}      
		for (i=0;i<=2*n-4;i++)
			LONGUEUR[i]=LONGUEUR[i]*pow(10.0,-POWER); 
	}              
}

/* Print tree edges wth their lengths into the Output file */
void PrintEdges (long int *ARETE, double *LONGUEUR, int EdgesNumber, int optionTR, double *CRITERION1, double *CRITERION2, int OptionFunction)
{  
	int i,j,NombreBlanc,a,b,c; 
	FILE *bootf;
	int addboot=0;
	
	if (optionTR!=2){ 
		fprintf(Output1,"\n\n              TREE EDGES WITH THEIR LENGTHS");
		if(p_boot == 0 && m_dataType == 2)
			fprintf(Output1," AND BOOTSTRAP SCORES (%d REPLICATES)",p_nbRep);
		if(p_boot == 1 && m_dataType == 2)
			fprintf(Output1," AND JACKKNIFE SCORES (%d REPLICATES)",p_nbRep);
		fprintf(Output1,"\n\n");	
	}
	else if (OptionFunction!=2) fprintf(Output1,"\n\n         RETICULOGRAM EDGES WITH THEIR LENGTHS  \tLSC      \tQ1\n\n");
	else fprintf(Output1,"\n\n         RETICULOGRAM EDGES WITH THEIR LENGTHS  \tLSC      \tQ2\n\n");
	
	fprintf(TreeFile,"\naretes = ");
	for (i=1;i<=EdgesNumber;i++){
		if(i!=1) fprintf(TreeFile,",");
		fprintf(TreeFile,"%d,%d",ARETE[2*i-1],ARETE[2*i-2]);
	}
	
	fprintf(TreeFile,"\nlongueur = ");
	for (i=1;i<=EdgesNumber;i++){
		if(i!=1) fprintf(TreeFile,",");
		fprintf(TreeFile,"%lf",LONGUEUR[i-1]);
	}
	
	if(p_boot<2){
		bootf = fopen(p_outboot,"r");
		if(bootf==NULL)
			addboot=0;
		else
			addboot=1;
	}
	for (i=1;i<=EdgesNumber;i++)
	{  
		NombreBlanc=0;
		if (ARETE[2*i-2]<10) NombreBlanc=3;
		else if (ARETE[2*i-2]<100) NombreBlanc=2;
		else if (ARETE[2*i-2]<1000) NombreBlanc=1;
		for (j=1;j<=NombreBlanc;j++)
			fprintf(Output1," ");                   
		fprintf(Output1,"            %d--%d",ARETE[2*i-2], ARETE[2*i-1]);
		
		NombreBlanc=0;
		if (ARETE[2*i-1]<10) NombreBlanc=3;
		else if (ARETE[2*i-1]<100) NombreBlanc=2;
		else if (ARETE[2*i-1]<1000) NombreBlanc=1;
		for (j=1;j<=NombreBlanc;j++)
			fprintf(Output1," ");  
		if ((optionTR!=2)||(i<2*n-3)){  
			fprintf(Output1,"            %f",  LONGUEUR[i-1]);
			if(addboot){
				fscanf(bootf,"%d %d %d",&a,&b,&c);
				fprintf(Output1," -- %d %%",c);
			}
			fprintf(Output1,"\n");
		}
		if ((optionTR==2)&&(i>=2*n-3)) fprintf(Output1,"            %f    \t%f    \t%f\n", LONGUEUR[i-1],CRITERION1[i-2*n+3],CRITERION2[i-2*n+3]);   
	}
	if ((optionTR==2)&&(EdgesNumber==2*n-3)) 
		fprintf(Output1,"\n\nNo extra edges have been added to the tree !");
	if (optionTR==2) 
		fprintf(Output1,"\n\nLeast-squares coefficient (LSC) and Q1 or Q2 values are given: first, for the starting additive tree, then, for each of extra edges added to it (if any)");
}


/* Reticulogram Reconstruction */
/* 
Input:
DISS - given dissimilarity matrix
D - additive distance matrix representing the starting tree
W - matrix of weights associated with dissimilarity entries (if no weights are considered, all Wij = 1, 0<i,j<=n)
OptionFunction - defines the stopping rule for the algorithm; possible values: 0 (a fixed number of edges should be placed into the additive tree), 1 (minimum of Q1 is reached), and 2 (minimum of Q2 is reached)
Iternumber - defines the maximum number of edges to be placed in the additive tree; Iternumber = n*(n-1)/2 if the value of OptionFunction is not 0.

Output:	 
ARETE - list of the reiculogram edges
LONGUEUR - list of the lengths of the reticulogram egges from ARETE
ReticulationsNumber - number of edges in the reticulogram
CRITERION1 - list of fit improvements according to the Q1 criterion for each edge (at each iteration) added to the starting additive tree
CRITERION2 - list of fit improvements according to the Q2 criterion for each edge (at each iteration) added to the starting additive tree  
*/ 
void RETICULATIONS (int n, double **DISS, double **D, long int *ARETE, double *LONGUEUR, int OptionFunction, double **W, int Iternumber, int *ReticulationsNumber, double *CRITERION1, double *CRITERION2, int unproc)
{
	int i,j,k,p,P,*X;
	double S,DIS,DIS1,*L,LON1,LON2,**DIST,FUNCTION_old,EQminGLold;
	int i1,j1,i2,j2,a,b,Aopt,Bopt,a10,**EDGES,NUM,NUMnumber,*NUM1,*A,iteration;
	double Sum,Sum1,Sum2,MIN,*A1,*LIMIT,EQ,EQmin,EQminGL,Lopt,LoptGL,Lopt1,FUNCTION,*NUMw,NUMnumberW;
	double *Tloptgl, *Teqmingl;; /* Tableau des EQminGL et LoptGL necessaire à la parallelisation */
	int *Taopt, *Tbopt; /* Tableau des aopt et bopt necessaire à la parallelisation */



	X=(int *)malloc((n+1)*sizeof(int));  
	L=(double *)malloc((n+1)*sizeof(double));
	double **Tree, **Path;
	Tree = (double**)malloc((2*n-2)*(sizeof(double*)));
	Path = (double**)malloc((2*n-2)*(sizeof(double*)));
	for (i=0;i<=2*n-3 ;i++ ){
		Path[i] = (double*)malloc( 3*(sizeof(double)));
		Tree[i] = (double*)malloc( 3*(sizeof(double))); 
	}
	DIST=(double **)malloc((2*n-1)*sizeof(double*));
	A=(int *)malloc((n*(n-1)+1)*sizeof(int)); 
	NUM1=(int *)malloc((n*(n+1)/2+1)*sizeof(int));
	NUMw=(double *)malloc((n*(n+1)/2+1)*sizeof(double));    
	EDGES=(int **)malloc((2*n-1)*sizeof(int*));  
	A1=(double *)malloc((n*(n+1)/2+1)*sizeof(double));
	LIMIT=(double *)malloc((n*(n+1)/2+1)*sizeof(double));

	for (i=0;i<=2*n-2;i++){
		DIST[i]=(double*)malloc((2*n-1)*sizeof(double));
		EDGES[i]=(int*)malloc((2*n-1)*sizeof(int));

		if ((DIST[i]==NULL)||(EDGES[i]==NULL)){ 
			printf("Data matrix is too large"); 
		exit(5);
		} 
	}
	/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1 || unproc==1){
		Taopt=(int *)malloc((numprocs+1)*sizeof(int)); 
		Tbopt=(int *)malloc((numprocs+1)*sizeof(int)); 
		Tloptgl=(double *)malloc((numprocs+1)*sizeof(double));
		Teqmingl=(double *)malloc((numprocs+1)*sizeof(double));
	}

	EQminGL=(double)LONG_MAX;

	i=1; j=n;
	odp(D,X,&i,&j);

	for (i=1;i<=n;i++){ 
		for (j=1;j<=n;j++) 
			DIST[i][j]=D[i][j];
	}

	/* Filling out the distance matrix DIST between vertices of the tree (then reticulogram) */
	L[1]=D[X[1]][X[2]];
	Path[1][0]=X[1];
	Path[1][1]=X[2];

	p=0;
	P=1;

	for(k=2;k<=n-1;k++){
		DIS=(D[X[1]][X[k]]+D[X[1]][X[k+1]]-D[X[k]][X[k+1]])/2;
		DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;
		if (DIS1<=0.00001) 
			DIS1=0.00001;
		S=0.0;
		i=0;

		if (DIS>0.00001){
			while (S<DIS-0.00001){
				i=i+1;
				S=S+L[i];
			}
		}
		else { 
			DIS=0; i=1; 
		}

		Tree[p+1][0]=n+k-1;
		Tree[p+1][1]=Path[i][1];
		Tree[p+1][2]=S-DIS;
		if (DIS<=0.00001) 
			Tree[p+1][2]=L[1];
		if (Tree[p+1][2]<0) 
			Tree[p+1][2]=0;

		for (j=i+1;j<=P;j++){
			Tree[p+j-i+1][0]=Path[j][0];
			Tree[p+j-i+1][1]=Path[j][1];
			Tree[p+j-i+1][2]=L[j];
			if (L[j]<0)
				L[j]=0; 
		}
		p=p+P-i+1;

		/* Filling out DIST in the iteration k*/
		LON1=DIS-(S-L[i]);
		LON2=S-DIS;
		if (DIS<=0.00001) { 
			LON1=0.0; LON2=L[1]; 
		} 

		for (j=1;j<=k;j++){    
			if (DIST[(int)Path[i][1]][X[j]]<DIST[(int)Path[i][0]][X[j]])      
				DIST[X[j]][n+k-1]=DIST[n+k-1][X[j]]=DIST[(int)Path[i][1]][X[j]]+LON2;
			else                      
				DIST[X[j]][n+k-1]=DIST[n+k-1][X[j]]=DIST[(int)Path[i][0]][X[j]]+LON1;
			if (j<=k-2){
				if (DIST[(int)Path[i][1]][n+j]<DIST[(int)Path[i][0]][n+j])      
					DIST[n+j][n+k-1]=DIST[n+k-1][n+j]=DIST[(int)Path[i][1]][n+j]+LON2;
				else                      
					DIST[n+j][n+k-1]=DIST[n+k-1][n+j]=DIST[(int)Path[i][0]][n+j]+LON1; 
			}     

		}
		DIST[n+k-1][n+k-1]=0.0;
		DIST[n+k-1][X[k+1]]=DIST[X[k+1]][n+k-1]=DIS1;
		for (j=n+1;j<=n+k-2;j++)
			DIST[j][X[k+1]]=DIST[X[k+1]][j]=DIST[n+k-1][j]+DIS1;

		/* End of the filling out DIST in the iteration k*/
		Path[i][1]=n+k-1;
		Path[i+1][0]=n+k-1;
		Path[i+1][1]=X[k+1];
		L[i]=L[i]+DIS-S;
		L[i+1]=DIS1;
		if (DIS<=0.00001) L[1]=0.00001; 
			P=i+1;
	}

	for (i=1;i<=P;i++) {
		Tree[p+i][0]=Path[i][0];
		Tree[p+i][1]=Path[i][1];
		Tree[p+i][2]=L[i];
	}

	for (i=1;i<=2*n-3;i++){
		if (fabs(Tree[i][2]-0.00001)<=0.00001)
			Tree[i][2]=0.0;
		ARETE[2*i-2]=Tree[i][0];
		ARETE[2*i-1]=Tree[i][1];
		LONGUEUR[i-1]=Tree[i][2];     
	} 

	EQ=0.0;
	for (i1=1;i1<=n-1;i1++){
		for (j1=i1+1;j1<=n;j1++)
			EQ=EQ+W[i1][j1]*(DISS[i1][j1]-DIST[i1][j1])*(DISS[i1][j1]-DIST[i1][j1]);
	}
	if (OptionFunction!=2)  
		FUNCTION=sqrt(EQ)/(n*(n-1)/2-2*n+3);
	else 
		FUNCTION=EQ/(n*(n-1)/2-2*n+3);
	EQminGLold=EQ; 
	CRITERION1[0]=EQ; 
	CRITERION2[0]=FUNCTION;
	for (i=1;i<=2*n-2;i++){
		for (j=1;j<=2*n-2;j++)
			EDGES[i][j]=0;
	}
	for (i=1;i<=2*n-3;i++)
		EDGES[ARETE[2*i-2]][ARETE[2*i-1]]=EDGES[ARETE[2*i-1]][ARETE[2*i-2]]=1;
	
	/* Main loop. Placement of new edges in the submitted additive tree.*/
	for (iteration=1;iteration<=Iternumber;iteration++){  
		if (n==2)
			break;
		a10=1;  
		/* Cas d'un processeur */
		if(numprocs ==1 || unproc==1){
			for (a=1;a<=2*n-3;a++){   	
				for (b=a+1;b<=2*n-2;b++){      
					if ((EDGES[a][b]==0)){       
						P=0; A1[0]=DIST[a][b]+1;
						for (i=1;i<=n-1;i++){ 
							for (j=i+1;j<=n;j++){
								if ((DIST[i][j]-DIST[a][i]-DIST[b][j])>(DIST[i][j]-DIST[a][j]-DIST[b][i])){
									MIN = DIST[i][j]-DIST[a][i]-DIST[b][j]; i2=i; j2=j; 
								}
								else{ 
									MIN = DIST[i][j]-DIST[a][j]-DIST[b][i]; i2=j; j2=i; 
								}     
								if (MIN>0.000001) {
									P=P+1; A[2*P-1]=i2; A[2*P]=j2; A1[P]=MIN; j1=P-1;
									while (A1[j1]<MIN){
										A[2*j1+1]=A[2*j1-1]; A[2*j1+2]=A[2*j1]; A1[j1+1]=A1[j1];
										A[2*j1-1]=i; A[2*j1]=j; A1[j1]=MIN; j1=j1-1;
									} 
								}
							} 
						}
						if (P>0) {
							NUM=1; NUM1[1]=1; NUMw[1]=W[A[1]][A[2]]; LIMIT[1]=A1[1]; 
							for (j1=2;j1<=P;j1++){ 
								if (A1[j1]<A1[j1-1]-0.000001){
									NUM=NUM+1; LIMIT[NUM]=A1[j1]; NUM1[NUM]=0; NUMw[NUM]=0;
								}
								NUM1[NUM]=NUM1[NUM]+1;
								NUMw[NUM]=NUMw[NUM]+W[A[2*j1-1]][A[2*j1]];
							} 
							LIMIT[NUM+1]=-0.000001; 
							Sum=Sum1=0.0; NUMnumber=0; NUMnumberW=0.0;
							for (i=1;i<=NUM;i++){
								NUMnumber=NUMnumber+NUM1[i]; NUMnumberW=NUMnumberW+NUMw[i];
								for (j=NUMnumber+1-NUM1[i];j<=NUMnumber;j++){
									i2=A[2*j-1]; j2=A[2*j];
									if (DIST[i2][a]+DIST[j2][b]-DIST[i2][b]-DIST[j2][a]>=0.0) { 
										j1=i2; i2=j2; j2=j1; 
									}
									Sum=Sum+W[i2][j2]*(DISS[i2][j2]-DIST[i2][a]-DIST[j2][b]);
								}
								{        
									Sum2=0.0; Sum1=0.0; double Sum0=0.0;         
									for (j=1; j<=NUMnumber; j++){
										i2=A[2*j-1]; j2=A[2*j];
										if (DIST[i2][a]+DIST[j2][b]-DIST[i2][b]-DIST[j2][a]>=0.0) {
											j1=i2; i2=j2; j2=j1; 
										}
										Sum2=Sum2+(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i])*(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i]);
										Sum1=Sum1+(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i+1]+0.000001)*(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i+1]+0.000001);         
										Sum0=Sum0+(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+Sum/NUMnumber)*(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+Sum/NUMnumber);
									}
								}
								if ((NUMnumberW!=0.0)&&((Sum/NUMnumberW>=LIMIT[i+1]+0.000001)&&(Sum/NUMnumberW<=LIMIT[i]))) { 
									Lopt1=Sum/NUMnumberW; 
								}
								else {
									if (Sum1<Sum2) { 
										Lopt1=LIMIT[i+1]+0.000001;
									}
									else { 
										Lopt1=LIMIT[i]; 
									}
								}
								EQ=0.0; i2=n*(n-1)/2; double EQ1=0.0;
								for (i1=1;i1<=n-1;i1++){
									for (j1=i1+1;j1<=n;j1++){         
										if (((DIST[i1][a]+DIST[j1][b]+Lopt1)<DIST[i1][j1]))
											EQ1=EQ1+W[i1][j1]*(DISS[i1][j1]-DIST[i1][a]-DIST[j1][b]-Lopt1)*(DISS[i1][j1]-DIST[i1][a]-DIST[j1][b]-Lopt1);
										else if ((DIST[i1][b]+DIST[j1][a]+Lopt1)<DIST[i1][j1]) 
											EQ1=EQ1+W[i1][j1]*(DISS[i1][j1]-DIST[i1][b]-DIST[j1][a]-Lopt1)*(DISS[i1][j1]-DIST[i1][b]-DIST[j1][a]-Lopt1); 
										else {
											EQ=EQ+W[i1][j1]*(DISS[i1][j1]-DIST[i1][j1])*(DISS[i1][j1]-DIST[i1][j1]); i2=i2-1;
										}
									}
								}
								EQ=EQ+EQ1;
								if ((i==1)||(EQmin>EQ)){ 
									Lopt=Lopt1; EQmin=EQ;
								}                 
							}  
							if (((a10==1)&&(P>0))||(EQminGL>EQmin)){
								EQminGL=EQmin;
								LoptGL=Lopt;
								Aopt=a;
								Bopt=b;
								a10=0;
							} 
						}   
					}
				}
			}
		}
		/* Cas de plusieurs processeurs */
		else{
			for (a=1;a<=2*n-3;a++){   
				if(a%numprocs==myid){
					for (b=a+1;b<=2*n-2;b++){      
						if ((EDGES[a][b]==0)){       
							P=0; A1[0]=DIST[a][b]+1;
							for (i=1;i<=n-1;i++){ 
								for (j=i+1;j<=n;j++){
									if ((DIST[i][j]-DIST[a][i]-DIST[b][j])>(DIST[i][j]-DIST[a][j]-DIST[b][i])){
										MIN = DIST[i][j]-DIST[a][i]-DIST[b][j]; i2=i; j2=j; 
									}
									else{ 
										MIN = DIST[i][j]-DIST[a][j]-DIST[b][i]; i2=j; j2=i; 
									}     
									if (MIN>0.000001) {
										P=P+1; A[2*P-1]=i2; A[2*P]=j2; A1[P]=MIN; j1=P-1;
										while (A1[j1]<MIN){
											A[2*j1+1]=A[2*j1-1]; A[2*j1+2]=A[2*j1]; A1[j1+1]=A1[j1];
											A[2*j1-1]=i; A[2*j1]=j; A1[j1]=MIN; j1=j1-1;
										} 
									}
								} 
							}
							if (P>0) {
								NUM=1; NUM1[1]=1; NUMw[1]=W[A[1]][A[2]]; LIMIT[1]=A1[1]; 
								for (j1=2;j1<=P;j1++){ 
									if (A1[j1]<A1[j1-1]-0.000001){
										NUM=NUM+1; LIMIT[NUM]=A1[j1]; NUM1[NUM]=0; NUMw[NUM]=0;
									}
									NUM1[NUM]=NUM1[NUM]+1;
									NUMw[NUM]=NUMw[NUM]+W[A[2*j1-1]][A[2*j1]];
								} 
								LIMIT[NUM+1]=-0.000001; 
								Sum=Sum1=0.0; NUMnumber=0; NUMnumberW=0.0;
								for (i=1;i<=NUM;i++){
									NUMnumber=NUMnumber+NUM1[i]; NUMnumberW=NUMnumberW+NUMw[i];
									for (j=NUMnumber+1-NUM1[i];j<=NUMnumber;j++){
										i2=A[2*j-1]; j2=A[2*j];
										if (DIST[i2][a]+DIST[j2][b]-DIST[i2][b]-DIST[j2][a]>=0.0) { 
											j1=i2; i2=j2; j2=j1; 
										}
										Sum=Sum+W[i2][j2]*(DISS[i2][j2]-DIST[i2][a]-DIST[j2][b]);
									}
									{        
										Sum2=0.0; Sum1=0.0; double Sum0=0.0;         
										for (j=1; j<=NUMnumber; j++){
											i2=A[2*j-1]; j2=A[2*j];
											if (DIST[i2][a]+DIST[j2][b]-DIST[i2][b]-DIST[j2][a]>=0.0) {
												j1=i2; i2=j2; j2=j1; 
											}
											Sum2=Sum2+(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i])*(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i]);
											Sum1=Sum1+(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i+1]+0.000001)*(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+LIMIT[i+1]+0.000001);         
											Sum0=Sum0+(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+Sum/NUMnumber)*(DIST[i2][a]+DIST[j2][b]-DISS[i2][j2]+Sum/NUMnumber);
										}
									}
									if ((NUMnumberW!=0.0)&&((Sum/NUMnumberW>=LIMIT[i+1]+0.000001)&&(Sum/NUMnumberW<=LIMIT[i]))) { 
										Lopt1=Sum/NUMnumberW; 
									}
									else {
										if (Sum1<Sum2) { 
											Lopt1=LIMIT[i+1]+0.000001;
										}
										else { 
											Lopt1=LIMIT[i]; 
										}
									}
									EQ=0.0; i2=n*(n-1)/2; double EQ1=0.0;
									for (i1=1;i1<=n-1;i1++){
										for (j1=i1+1;j1<=n;j1++){         
											if (((DIST[i1][a]+DIST[j1][b]+Lopt1)<DIST[i1][j1]))
												EQ1=EQ1+W[i1][j1]*(DISS[i1][j1]-DIST[i1][a]-DIST[j1][b]-Lopt1)*(DISS[i1][j1]-DIST[i1][a]-DIST[j1][b]-Lopt1);
											else if ((DIST[i1][b]+DIST[j1][a]+Lopt1)<DIST[i1][j1]) 
												EQ1=EQ1+W[i1][j1]*(DISS[i1][j1]-DIST[i1][b]-DIST[j1][a]-Lopt1)*(DISS[i1][j1]-DIST[i1][b]-DIST[j1][a]-Lopt1); 
											else {
												EQ=EQ+W[i1][j1]*(DISS[i1][j1]-DIST[i1][j1])*(DISS[i1][j1]-DIST[i1][j1]); i2=i2-1;
											}
										}
									}
									EQ=EQ+EQ1;
									if ((i==1)||(EQmin>EQ)){ 
										Lopt=Lopt1; EQmin=EQ;
									}                 
								}  
								if (((a10==1)&&(P>0))||(EQminGL>EQmin)){
									Teqmingl[myid]=EQmin;Tloptgl[myid]=Lopt;Taopt[myid]=a;Tbopt[myid]=b;EQminGL=EQmin;LoptGL=Lopt;Aopt=a;Bopt=b;a10=0;
								} 
							}   
						}
					}
				}
			
			}
			/* Synchronisation des resultats */
			if(myid>0){
				MPI_Send( &LoptGL, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
				MPI_Send( &EQminGL, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
				MPI_Send( &Aopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
				MPI_Send( &Bopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
				MPI_Recv( &LoptGL, 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &EQminGL, 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Aopt, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Bopt, 1, MPI_INT, myid-1,0, MPI_COMM_WORLD, &status );
				if(myid<numprocs-1){
					MPI_Send( &LoptGL, 1, MPI_DOUBLE,myid+1, 0, MPI_COMM_WORLD );
					MPI_Send( &EQminGL, 1, MPI_DOUBLE,myid+1, 0, MPI_COMM_WORLD );
					MPI_Send( &Aopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
					MPI_Send( &Bopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				}
			}
			else{
				for(i=1; i<numprocs;i++){
					MPI_Recv( &Tloptgl[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Teqmingl[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Taopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &Tbopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
				}
				for(i=1; i<numprocs;i++){
					if( Teqmingl[i]< EQminGL){
						Aopt = Taopt[i];Bopt = Tbopt[i];EQminGL = Teqmingl[i];LoptGL = Tloptgl[i];
					}
					if((Teqmingl[i]==EQminGL) && (Taopt[i] < Aopt )){
						Aopt = Taopt[i];Bopt = Tbopt[i];EQminGL = Teqmingl[i];LoptGL = Tloptgl[i];
					}	
				}
				MPI_Send( &LoptGL, 1, MPI_DOUBLE,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &EQminGL, 1, MPI_DOUBLE,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &Aopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &Bopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			}
		}
		if (a10==1)
			break;
		else {
			FUNCTION_old=FUNCTION; 
			if ((OptionFunction!=2)&&(iteration<n*(n-1)/2-2*n+3))
				FUNCTION=sqrt(EQminGL)/(n*(n-1)/2-2*n+3-iteration);
			else if (iteration<n*(n-1)/2-2*n+3) 
				FUNCTION=EQminGL/(n*(n-1)/2-2*n+3-iteration);
			if ((Iternumber>=(2*n-2)*(2*n-2-1)/2-2*n+3)&&(iteration==n*(n-1)/2-2*n+3)) 
				break;
			if ((Iternumber>=(2*n-2)*(2*n-2-1)/2-2*n+3)&&(FUNCTION_old<=FUNCTION))
				break;   
			if ((Iternumber==(2*n-2)*(2*n-2-1)/2-2*n+3)&&(fabs(EQminGL-EQminGLold)<=0.00000001))
				break;
			EQminGLold=EQminGL;
			CRITERION1[iteration]=EQminGL;
			CRITERION2[iteration]=FUNCTION;
			ARETE[4*n-6+2*iteration-2]=Aopt;
			ARETE[4*n-6+2*iteration-1]=Bopt;
			LONGUEUR[2*n-4+iteration]=LoptGL;
			EDGES[Aopt][Bopt]=EDGES[Bopt][Aopt]=1;

			/* Updating the matrix DIST after addition of a new edge*/

			for (i=1;i<=2*n-3;i++){ 
				for (j=i+1;j<=2*n-2;j++){
					{
						if (DIST[i][Aopt]+DIST[j][Bopt]+LoptGL<DIST[i][j]) 
							DIST[i][j]=DIST[j][i]=DIST[i][Aopt]+DIST[j][Bopt]+LoptGL;
						else if (DIST[i][Bopt]+DIST[j][Aopt]+LoptGL<DIST[i][j]) 
							DIST[i][j]=DIST[j][i]=DIST[i][Bopt]+DIST[j][Aopt]+LoptGL; 
					}
				}              
			}   
		}
		if ((2*n-2)*(2*n-2-1)/2-2*n+3-iteration==0)   { 
			*ReticulationsNumber=(2*n-2)*(2*n-2-1)/2; break;  
		}

	}

	if (*ReticulationsNumber==0) 
		*ReticulationsNumber=2*n-3+(iteration-1);

	for (i=1;i<=n;i++){ 
		for (j=1;j<=n;j++)
			D[i][j]=DIST[i][j];
	}
	/*if (myid == 0) {
		endwtime = MPI_Wtime();
		totaltime = totaltime + endwtime - startwtime;
	}*/

	free(X);
	free(Tree);
	free(L);
	free(Path);

	free(A);
	free(NUM1);
	free(NUMw);
	free(A1);
	free(LIMIT);

	for (i=0;i<=2*n-2;i++){
		free(DIST[i]); free(EDGES[i]);
	}
	free(DIST); free(EDGES); 
}



/* Tree reconstruction from incomplete matrices*/


/* Main function calling the tree reconstruction method from incomplete matrices*/
void GL_Main(double **DD, int nbObj, int methode, double **distanceArbre, double *RESULTATS, long int *ARETE, double *LONGUEUR, double **W)
{
	int add, complete;	   	  
	int i, j, T, nbmiss, nbultra, nbcycles, sum, sum2=0, **B2, *aa; 
	double  *max, **m;     	
	int tt, k; 
	int NbInc=0;
	double **DI;
	int error, Iternum=50;
	int FLAG1;
	
	n=nbObj;
	if (methode!=4) {
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{
				k=1+floor1((n-0.5*i)*(i-1)+j-i-1);    
				if ((DD[i][j]<0.0)&&(DD[i][j]!=-99)) 
				{printf("Dissimilarity matrix contains negative values not equal to -99 (value marking missing data). Computation is not possible for this method."); 
				exit(5); }	  			
			}
		}	
	}
	
	DI=(double **) malloc((n+1)*sizeof(double*)); 
	for (i=0;i<=n;i++)
	{
		DI[i]=(double*)malloc((n+1)*sizeof(double));
		if (DI[i]==NULL)
		{
			printf("Data matrix is too large"); 
			exit(5); 
		}
	}
	
	{ 		
		Dist=(double **) malloc((n+1)*sizeof(double*)); 
		for (i=0;i<=n;i++)
		{
			Dist[i]=(double*)malloc((n+1)*sizeof(double));
			if (Dist[i]==NULL)
			{
				printf("Data matrix is too large"); exit(5); 
			}
		}	
	}
	
	for (i=1;i<=n;i++)
	{
		for (j=i+1;j<=n;j++)
		{
			k=1+floor1((n-0.5*i)*(i-1)+j-i-1);    			
			DI[i][j]=DD[i][j];
			
			DI[j][i]=DI[i][j];
			if (methode==1) { Dist[j-1][i-1]=Dist[i-1][j-1]=DI[j][i]; }
			
			if (Dist[i-1][j-1]<0)	NbInc++;
			else if (Dist[i-1][j-1]>Dmax) Dmax=Dist[i-1][j-1];
		}
		if (methode==1) { DI[i][i]=Dist[i-1][i-1]=0.0; }
		else DI[i][i]=0.0;
	}
	
	/* Test for avoiding entire rows and columns with missing entries*/
	FLAG1;
	for (i=1;i<=n;i++)
	{
		FLAG1=1;
		for (j=1;j<=n;j++)
		{			   
			if ((i!=j)&&(DI[i][j]!=-99)) { FLAG1=0; break; }
		}
		if (FLAG1==1) { printf("Too many missing values. You have to have at least one value by row and by colomn.");
		exit(5); } 
	}


	/*Selection of a tree reconstruction method*/
	if (methode==1) /*Triangles method*/
	{
		N=nbObj;
	
		error=PSM(Dist,N);

		if (error==-1) {printf("Too many missing entrees for the Triangles method"); 
		exit(5);} 
		
		for (i=1;i<=2*N-3;i++)
		{	
			ARETE[2*i-2]=i;
			ARETE[2*i-1]=Tree[i-1]+1;
			LONGUEUR[i-1]=Long[i-1];  
		}
		
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{
				k=1+floor1((n-0.5*i)*(i-1)+j-i-1); 
				distanceArbre[i][j]=distanceArbre[j][i]=Dist[i-1][j-1];
			}
			distanceArbre[i][i]=0.0;	     
		}		
		
		RESULTATS[0]=RESULTATS[1]=RESULTATS[2]=RESULTATS[3]=0.0;
		for (i=1;i<=n-1;i++)
		{
			for (j=i+1;j<=n;j++)
			{  
				double s=fabs(DI[i][j]-Dist[i-1][j-1]);    
				if (DI[i][j]!=-99) RESULTATS[0]=RESULTATS[0]+s*s;        
				if (DI[i][j]!=-99) RESULTATS[1]=RESULTATS[1]+s;
				if ((DI[i][j]!=-99)&&(s>RESULTATS[2])) RESULTATS[2]=s;
			}
		} 
		RESULTATS[1]=RESULTATS[1]/(n*(n-1)/2.0);
		for (i=1;i<=2*n-3;i++)
			RESULTATS[3]=RESULTATS[3]+LONGUEUR[i-1];	
	}
	
	
	else if ((methode==2)||(methode==3)) /*Ultrametric and Additive procedures*/
	{
		T=n;	
		max = (double*)malloc((T+1) * sizeof(double));
		m = (double**)malloc((T+1) * sizeof(double *));
		B2 = (int**)malloc((T+1) * sizeof(int *));		
		if ((max== NULL) ||(m == NULL) || (B2 == NULL))
		{ printf("Data matrix is too large"); 
		exit(5);} 	
		
		for (i=0;i<T+1;i++)
		{ 			
			m[i] = (double*)malloc((T+1) * sizeof(double));
			B2[i] = (int*)malloc((T+1) * sizeof(int));
			if ((m[i] == NULL) || (B2[i] == NULL))
			{ printf("Data matrix is too large"); 
			exit(5);} 
		}				
		
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{
				k=1+floor1((n-0.5*i)*(i-1)+j-i-1);    
				m[i][j]=m[j][i]=DD[i][j];
			}			
			m[i][i]=0.0;
		}
		
		nbmiss=0;
		complete=1;	
		tt=(T * (T - 1))/2;
		for (i=1; i<=T; i++)
		{
			for (j=1; j<=T; j++)
			{if (m[i][j]!=miss) 
				B2[i][j]=1;
				else 
			{ B2[i][j]=0; nbmiss++;}}
		}
		if (nbmiss>0) complete=0;    
		
		aa=(int*)malloc((T+1) * sizeof(int));
		
		if (methode==2) add=0; 
		if (methode==3) add=1;  
		
		if (nbmiss > 0)
		{
			nbultra=0;
			nbcycles=1;
			if (add) Additif(&sum, &sum2, T, B2, m, max);
			else Ultra2(&sum, &sum2, T, B2, m, max);
		}	 
		while (sum2>0)
		{
			nbcycles=nbcycles+1;
			if (!add) Ultra2(&sum, &sum2, T, B2, m, max);
			if (add)
			{
				if (sum2==sum)
				{
					nbultra=nbultra+1;
					alea(T, aa);
					Ultra1(&sum, &sum2, T, B2, m, max, aa);
				}
				Additif(&sum, &sum2, T, B2, m, max);
			}
		}
		
		
		/*method MW to complete Additive and Ultrametric*/
		for (i=1; i<=n; i++)
		{
			for (j=i+1; j<=n; j++)
			{ 
				W[i][j]=W[j][i]=1.0;
			}   
			W[i][i]=1.0;
		} 
		W[0][0]=1.0;
		
		for (i=1; i<=n; i++)
		{
			for (j=1; j<=n; j++)
				if (m[i][j] <= 0.0) m[i][j] = 0.0005;
			
		}
		
		error=parcour211(m,W,Dist,&Iternum,ARETE,LONGUEUR); 
		if (error==-1) {printf("encore un probleme...");exit(5);}
		
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{
				k=1+floor1((n-0.5*i)*(i-1)+j-i-1); 
				distanceArbre[i][j]=distanceArbre[j][i]=Dist[i][j];
			}
			distanceArbre[i][i]=0.0;	     
		}		
		compute_criteres11(DI, Dist, RESULTATS, LONGUEUR, n);
	}
	
	else if ((methode==4)||(methode==5)) /* MW-modified method or MW* */
	{
		
		int FLAG=0;
		for (i=1; i<=n; i++)
		{
  			for (j=i+1; j<=n; j++)
			{  		  		
  				if ((DI[i][j]==-99.00)&&(DI[j][i]==-99.00)) {W[i][j]=W[j][i]=0.0;FLAG=1;}
  				else if ((DI[i][j]==-99.00)&&(DI[j][i]!=-99.00)) {DI[i][j]=DI[j][i]; W[i][j]=W[j][i]=1.0;}
  				else if ((DI[i][j]!=-99.00)&&(DI[j][i]==-99.00)) {DI[j][i]=DI[i][j]; W[i][j]=W[j][i]=1.0;}
  				else {W[i][j]=W[j][i]=1.0;}
			}   
			DI[i][i]=0.0;
			W[i][i]=1.0;
		}
		DI[0][0]=0; W[0][0]=0;
		
		if ((FLAG==1)&&(methode==5)) FLAG=tryout(DI, W);     
		
		if (n==2)
		{
			Dist[1][2]=DI[1][2];
			Dist[2][1]=DI[1][2];
			Dist[1][1]=0.0;
			Dist[2][2]=0.0;  
		}      
		else 
		{    
			Iternum=50; 
			DI[0][0]=1;
			
			printf("\nAvant traitement principal ");
			if (FLAG==0) error=parcour211(DI,W,Dist,&Iternum,ARETE,LONGUEUR); 						   
			else error=parcour2MisVal(DI,W,Dist,&Iternum,ARETE,LONGUEUR);
			
			if (error==-1) { printf("j'te jure!!!"); exit(5);}
		}   
		
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{
				k=1+floor1((n-0.5*i)*(i-1)+j-i-1); 
				distanceArbre[i][j]=distanceArbre[j][i]=Dist[i][j];
			}	     
		}		
		compute_criteres11(DD, Dist, RESULTATS, LONGUEUR, n);
	}
	
	

	/* memory liberation*/
	for (i=0;i<=n;i++)
	{ 
		free(DI[i]);
		free(Dist[i]);
		if ((methode==2)||(methode==3)) free(m[i]);
		if ((methode==2)||(methode==3)) free(B2[i]);
	}    
	free(DI);
	free(Dist);
	if ((methode==2)||(methode==3)) free(m);
	if ((methode==2)||(methode==3)) free(B2);
	if ((methode==2)||(methode==3)) free(max);
	return;				
}

/* Triangles method functions*/
void EditDist(double **D, int N)
{	int i,j;
	
	printf("Distance %s\n",Nom);
  	printf("N = %3d\n",N);
	for (i=0;i<N;i++)
	{ 	
		printf("%3d ",i+1);
		for (j=0;j<N;j++)
			printf("%6.2f",D[i][j]);
		printf("\n");
	}	
	printf("\n"); return;
}

void Chemin(int *T, int t, int u)
{	int k;
	
	for (k=0;k<2*N-2;k++) Ch[k]=0;
	k=t;
	while (k<2*N-3) { Ch[k]=1; k=T[k]; }
	k=u;
	while (k<2*N-3) { Ch[k]+=1; k=T[k]; }
	return;
}

void Arete(double *path, int *x, int *y)
{	double Lon;
	int Ix,Iy;
	
	Lon=*path; Ix=*x; Iy=*y;
	while ( Long[Ix]+eps < Lon )
	{	Lon -= Long[Ix]; Ix=Tree[Ix];
		if (Ix==0 || Ch[Ix]!=1) return;
	}
	*x=Ix; *y=Tree[Ix]; *path=Lon; return;
}

int GrowTree(int s, int t, int u, int NbPl)
/* Expending the tree with adjusted length for the new edge*/
{	double Longt, Longu, Lon;
	int i, ii, j, jj, Ix, Iy, cx, cy, NbP, flag;
	
	Ns--; Tree[s]=Ns;
	Long[s] = (Dist[s][t]+Dist[s][u]-Dist[t][u])/2;
	if (Long[s]<0) Long[s]=0;
	Longt = (Dist[s][t]+Dist[t][u]-Dist[s][u])/2;
	Longu = (Dist[s][u]+Dist[t][u]-Dist[t][s])/2;
	if (Longt<0) { Longu += Longt; Longt=0; }
	if (Longu<0) { Longt += Longu; Longu=0; }
	Chemin(Tree,t,u);
	
	Lon=Longt; Ix=t; Iy=0; Arete(&Lon, &Ix, &Iy);
	if (Iy==0) 	
		{	Lon=Longu; Ix=u; Iy=0; Arete(&Lon, &Ix, &Iy); }
	if (Iy==0) { 			
		{ printf("Too many missing values in the data matrix"); 
		exit(5);} 		
	}
	
	Tree[Ns]=Iy; Long[Ns]=Long[Ix]-Lon;
	Tree[Ix]=Ns; Long[Ix]=Lon;
	
	cx=0; cy=0;
	for (i=0;i<NbPl;i++)
	{	ii=Pl[i]; flag=1;
		while (flag)
		{	if (ii==Ix) {Clx[cx]=Pl[i]; cx++; flag=0; }
			if (ii==Iy || ii==2*N-3) {Cly[cy]=Pl[i]; cy++; flag=0; }
			ii=Tree[ii];
		}
	}	Lon=0; NbP=0;
	
	for (i=0;i<cx;i++)
	{	ii=Clx[i]; if (Dist[s][ii]<0) continue;
		for (j=0;j<cy;j++)
		{	jj=Cly[j]; if (Dist[s][jj]<0) continue;
			Lon += Dist[s][ii]+Dist[s][jj]-Dist[ii][jj];
			NbP++;
	}	}
	Long[s]=Lon/NbP/2; if (Long[s]<0) Long[s]=0.;
	return 1;
}


void MajDist(int s, int NbPl)
/* Compute the distances in the tree between s and all the placed objects */
{	int i, ii, k;
	double Lon;
	
	for (i=0;i<NbPl;i++)
	{	ii=Pl[i]; Lon=0.;
		Chemin(Tree,s,ii);
		for (k=0;k<2*N-3;k++)
		{	if (Ch[k]==1) Lon += Long[k];
			Dist[s][ii]=Lon; Dist[ii][s]=Lon;
			
	}	}	return;
}

void IniTri(double **D, int *ps, int *pt, int *pu)
/* Initialization with a triangle metric of a minimum perimeter */
{	int i, j, k;
	double Lon, Peri;
	
	Peri=3*Dmax+1;
	for (i=0;i<N-2;i++)
		for (j=i+1;j<N-1;j++)
	{	if (D[i][j]<0) continue;
		for (k=j+1;k<N;k++)
		{	if (D[i][k]<0 || D[j][k]<0) continue;
			if (D[i][j]>D[i][k]+D[j][k] || D[i][j]<fabs(D[i][k]-D[j][k])) continue;
			Lon=D[i][j]+D[i][k]+D[j][k];
			if (Lon<Peri)
			{	Peri=Lon;
				*ps=i; *pt=j; *pu=k;
	}	}	}
	if (Peri<3*Dmax+1) return;
	{ printf("There is no metric triangle. Try another reconstruction method."); 
	exit(5);} 
	
	*ps=0; *pt=1; *pu=2;
	return;
}

int PSM(double **D, int N)
{	int i,j,k, ii,jj, s=0,t=1,u=2, *Mark, NbPl;
	double Lon, Slong=0., Lk;
	char prog[5]="PSM.";
	int error;
	
	Mark = (int*) malloc((N) * sizeof(int));
	Pl = (int*) malloc((N) * sizeof(int));				/* list of placed vertices */
	Long = (double*) malloc((2*N) * sizeof(double));		/* edge lengths  */
	Tree = (int*) malloc((2*N) * sizeof(int));			/* Edge i -- Tree[i]  */
	Ch = (int*) malloc((2*N) * sizeof(int));			/* Path index for the edges */
	Clx = (int*) malloc((N) * sizeof(int));			/* Leaves located on the one */
	Cly = (int*) malloc((N) * sizeof(int));			/* and on the other side of an internal vertex */
	for (i=0;i<N;i++)
		Mark[i]=0;
	
	/* Initialization with a metric triangle of a minimum perimeter */
	IniTri(D,&s,&t,&u);
	NbPl=3; Ns=2*N-3;
	Mark[s]=1; Mark[t]=1; Mark[u]=1;
	Pl[0]=s; Pl[1]=t; Pl[2]=u;
	Long[s] = (D[s][t]+D[s][u]-D[t][u])/2;
	Long[t] = (D[s][t]+D[t][u]-D[s][u])/2;
	Long[u] = (D[u][t]+D[s][u]-D[t][s])/2;
	Tree[s]=Ns; Tree[t]=Ns; Tree[u]=Ns;

	/* Adding a vertex nearest to a path */
	while (NbPl<N)
	{	
		Lon=2*Dmax; s=-1;
		for (k=0;k<N;k++){	
			if (Mark[k]==1) continue; /* already placed */
			for (i=0;i<NbPl-1;i++)
			{	ii=Pl[i]; if (D[k][ii]<0) continue;
				for (j=i+1;j<NbPl;j++)
				{	jj=Pl[j]; if (D[k][jj]<0) continue;
					Lk=D[k][ii]+D[k][jj]-D[ii][jj];
					if (Lk<Lon) { Lon=Lk; s=k; t=ii; u=jj; }
				}	
			}	
		}

		if (s<0) { 
			free(Mark); free(Pl); free(Long); free(Tree); free(Ch); free(Clx); free(Cly);
			printf("Too many missing values in the data matrix"); 
		exit(5);} 
		error=GrowTree(s,t,u,NbPl);
		
		if (error==-1) { 
			free(Mark); free(Pl); free(Long); free(Tree); free(Ch); free(Clx); free(Cly);
			printf("To many missing values in the data matrix"); 
		exit(5);} 
		MajDist(s,NbPl);
		Mark[s]=1; Pl[NbPl]=s; NbPl++;
	}
	
	free(Mark);
	free(Pl); 
	free(Ch);
	free(Clx);
	free(Cly);
	return 1;
}


/* Ultrametric and Additive Procedures*/

void LecDist(double **m, int **B2, double *max, int *nbmiss, int *complete, int *T)
{	
	int i,j,tt; 
	char **Et, Nom[20];
	int NbInc=0;
	int T_temp;
	FILE *FichDist;
	
	printf("Distance file ");
	gets(Nom);
	printf("%s\n",Nom);
	FichDist = fopen(Nom,"r"); assert(FichDist != NULL);
	
	fscanf(FichDist,"%d",&T_temp);	
	*T = T_temp;
	
	max = (double*)malloc((T_temp+1) * sizeof(double));
	m = (double**)malloc((T_temp+1) * sizeof(double *));
	B2 = (int**)malloc((T_temp+1) * sizeof(int *));		
	Et = (char**)malloc((T_temp+1) * sizeof(char *));
	assert(max != NULL && m != NULL && B2 != NULL && Et != NULL);
	
	for (i=0;i<T_temp+1;i++)
	{ 			
		m[i] = (double*)malloc((T_temp+1) * sizeof(double));
		B2[i] = (int*)malloc((T_temp+1) * sizeof(int));
		Et[i] = (char*)malloc(20);
		assert(m[i] != NULL && B2[i] != NULL && Et[i] != NULL);
	}
	
	for (i=1;i<T_temp+1;i++)
	{ 
		fscanf(FichDist,"%s",Et[i]);
		for (j=1;j<T_temp+1;j++)
		{ double cc;
			fscanf(FichDist,"%f",&cc);
		m[i][j]=cc;}
	}	
	
  	fclose(FichDist);	   
    	
	*nbmiss=0;
	*complete=1;	
	tt=(T_temp * (T_temp - 1))/2;
	for (i=1; i<=*T; i++)
	{
		for (j=1; j<=*T; j++)
		{if (m[i][j]!=miss) 
			B2[i][j]=1;
			else 
		{ B2[i][j]=0; *nbmiss++;}}
	}
	if (*nbmiss>0) *complete=0;    
} 

void alea (int T, int *aa)
{  
	int i, j, k;
	double x, *P;
	double XY=0.98674532;
	
	P = (double*)malloc((T+1) * sizeof(double));
	assert(P!=NULL);
	
	for (i=1;i<=T;i++)
	{
		double D2P31M=2147483647.0, D2P31=2147483648.0;                                                      
		double a;                                                                   
		a=16807 * XY;                           
		XY=a - floor1(a/D2P31M) * D2P31M;
		P[i]=XY/D2P31; 
	}
	
	x=1;
	for (i=1;i<=T;i++)
	{
		for (j=1;j<=T;j++)
		{
			if (P[j] < x) 
				{ x = P[j]; aa[i] = j; }		        		     
		}       
		for (k=1;k<=T;k++)
		{
			if (P[k] == x) P[k] = 1;
		}
		x=1;
	} 
}


void Ultra1(int *sum, int *sum2, int T, int **B2, double **m, double *max, int *aa)
{			
	int i, j, k, a;
	double min, val1, val2;
	int done;
	
	done = 0;
	*sum = 0;
	*sum2 = 0;
	min = 1000;
	for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[aa[i]][aa[j]]==0) *sum=*sum+1;
		} 
	}   
	for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[aa[i]][aa[j]]==1) goto label2;
			for (k=1;k<=T;k++)
                        {
				if ((k==i)||(k==j)) goto label1;									
				if ((B2[aa[i]][aa[k]]==0)||(B2[aa[j]][aa[k]]==0)) goto label1;																	
				val1=(m[aa[i]][aa[k]] + m[aa[k]][aa[i]])/2.0;
				val2=(m[aa[j]][aa[k]] + m[aa[k]][aa[j]])/2.0;
				if (val1>val2) max[aa[k]]=val1;          			   
                                else max[aa[k]]=val2;
				
				if (max[aa[k]] < min)
				{
					min=max[aa[k]];
					m[aa[i]][aa[j]]=min;
					m[aa[j]][aa[i]]=min;     
					done=1;
				}
				label1:                      a=1;
			}
			min=1000;
			if (done) goto label3;
			label2:                  a=1;
		}			
	}
        
	label3:     a=1;
	for (i=1;i<=T;i++)
	{
		for (j=1;j<=T;j++)
		{
			if (m[aa[i]][aa[j]]!=miss) B2[aa[i]][aa[j]]=1;			
		}
	}  			
	
        for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[aa[i]][aa[j]]==0) *sum2=*sum2+1;
		}
	} 		
} 

void Ultra2(int *sum, int *sum2, int T, int **B2, double **m, double *max)
{			
	int i, j, k, a;
	double min, val1, val2;
	
	*sum = 0;
	*sum2 = 0;
	min = 1000;
	for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[i][j]==0) *sum=*sum+1;
		} 
	}   
	for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[i][j]==1) goto label2;
			for (k=1;k<=T;k++)
                        {
				if ((k==i)||(k==j)) goto label1;									
				if ((B2[i][k]==0)||(B2[j][k]==0)) goto label1;																	
				val1=(m[i][k] + m[k][i])/2.0;
				val2=(m[j][k] + m[k][j])/2.0;
				if (val1>val2) max[k]=val1;          			   
                                else max[k]=val2;
				
				if (max[k] < min)
				{
					min=max[k];
					m[i][j]=min;
					
				}
				label1:                      a=1;
			}
			min=1000;
			label2:                  a=1;
		}
	}
        
	for (i=1;i<=T;i++)
	{
		for (j=1;j<=T;j++)
		{
			if (m[i][j]!=miss) B2[i][j]=1;			
		}
	}  			
	
        for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[i][j]==0) *sum2=*sum2+1;
		}
	} 		
} 

void Additif(int *sum, int *sum2, int T, int **B2, double **m, double *max)
{	
	int i, j, k, l, a;
	double min, val1, val2, val3, newval;
	
	*sum=0;
	*sum2=0;
	min=1000;
	for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[i][j]==0) *sum=*sum+1;
		} 
	}   
	
	for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[i][j]==1) goto label3;
			for (k=1;k<=T;k++)
			{									
				if ((k==i)||(k==j)) goto label2;
				if ((B2[i][k]==0)||(B2[j][k]==0)) goto label2;
				for (l=1;l<=T;l++)
				{
					if ((l==i)||(l==j)||(l==k)) goto label1;
					if ((B2[i][l]==0)||(B2[j][l]==0)||(B2[k,l]==0)) goto label1;
					val1=(m[k][l]+m[l][k])/2.0;
					val2=(m[i][l]+m[l][i])/2.0 + (m[j][k]+m[k][j])/2.0;
					val3=(m[i][k]+m[k][i])/2.0 + (m[j][l]+m[l][j])/2.0;
					if (val2>val3) max[l]=val2;
					else max[l]=val3;
					newval=max[l]-val1;											
					if (newval<min)
					{
						min=newval;
						m[i][j]=min;
					}
					label1:                                  a=1;
				}
				label2:                          a=1;
			}
			min=1000;
			label3:                  a=1;
		}
	}				
	
	for (i=1;i<=T;i++)
	{
		for (j=1;j<=T;j++)
		{
			if (m[i][j]!=miss) B2[i][j]=1;			
		}
	}  							     
        for (i=1;i<=T;i++)
        {
		for (j=1;j<=T;j++)
		{
			if (B2[i][j]==0) *sum2=*sum2+1;
		}
	} 										
}

/*Computing the values of criteria*/
void compute_criteres11(double **D,double **DA, double *RESULTATS, double *LONGUEUR, int nn)
{
	int i1,j1;
	double s, R[4];
	
	
	R[0]=0.0;
	R[1]=0.0;
	R[3]=0.0;
	R[2]=0.0;
	for (i1=1;i1<=nn-1;i1++)
	{
		for (j1=i1+1;j1<=nn;j1++)
		{   
			s=fabs(D[i1][j1]-DA[i1][j1]);    
			if (D[i1][j1]!=-99.0)
			{ R[0]=R[0]+s*s;        
				R[1]=R[1]+s;
			if (s>R[2]) R[2]=s; }
		}
	} 
	
	for (i1=1;i1<=2*nn-3;i1++)
		R[3]=R[3]+LONGUEUR[i1-1];
	
	
	RESULTATS[0]=R[0];
	RESULTATS[1]=R[1]/(n*(n-1)/2.0);
	RESULTATS[2]=R[2];
	RESULTATS[3]=R[3];
	
}

/*MW-modified method */
int parcour2MisVal(double **DISS,double **W,double **TM,int *Iternum,long int *ARETE, double *LONGUEUR)
{
	int a,i,j,i1,j1,iopt,jopt,*X,Iternumber,Option,k,k1,*Flag;
	double **TM1, **TMopt;
	double EQ,EQmin,MAX,*SumW,*Sum,MAX1;
	int error1;
	int leproc;/* Numero du processeur ayant les meilleurs resultats */
	int *Tiopt, *Tjopt,*Ta; /* Tableau des iopt, jopt et a necessaire à la parallelisation */
	double *Teqmin;	 /* Tableau des EQmin necessaire à la parallelisation */

	/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1){
		Tiopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Ta=(int *) malloc((numprocs+1)*sizeof(int)); 
		Tjopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Teqmin=(double*)malloc((numprocs+1)*sizeof(double)); 
	}

	EQmin=(double)LONG_MAX;

	X=(int *) malloc((n+1)*sizeof(int)); 
	Flag=(int *) malloc((n+1)*sizeof(int)); 
	TM1=(double **) malloc((n+1)*sizeof(double*));
	TMopt=(double **) malloc((n+1)*sizeof(double*));
	SumW=(double *) malloc((n+1)*sizeof(double));
	Sum=(double *) malloc((n+1)*sizeof(double));
	
	for (i=0;i<=n;i++)
	{
		TM1[i]=(double*)malloc((n+1)*sizeof(double));
		TMopt[i]=(double*)malloc((n+1)*sizeof(double));  
		if ((TM1[i]==NULL)||(TMopt[i]==NULL))
		{
			printf("Data matrix is too large"); 
			exit(5); 
		}  
	}    
	
	Iternumber=*Iternum;
	X[0]=1;
	a=0;
	
	for (i=1;i<=n;i++)
	{
		Sum[i]=0.0;
		for (j=1;j<=n;j++)
			Sum[i]=Sum[i]+W[i][j];
	}
	
	Option=1;
	/* Cas d'un processeur */
	if(numprocs==1){
		for (i=1;i<=n;i++){ 
			for (j=i+1;j<=n;j++)   
			{
				if (W[i][j]!=0.0)
				{ 
					X[1]=i;
					X[2]=j;                    
					
					for (k=1;k<=n;k++)
						Flag[k]=0;
					Flag[i]=Flag[j]=1;
					for (k=1;k<=n;k++)
					{ 
						if (Flag[k]!=1)   
						{
							SumW[k]=W[k][i];
						}
					}
					
					for (k1=2;k1<=n-1;k1++)
					{
						MAX=MAX1=0.0;
						for (k=1;k<=n;k++)
						{ 
							if (Flag[k]!=1)   
							{             
								SumW[k]=SumW[k]+W[k][X[k1]];
								if ((MAX<SumW[k])||((MAX==SumW[k])&&(MAX1<Sum[k]))){ MAX=SumW[k]; X[k1+1]=k; MAX1=Sum[k];}    
							}              
						}
						Flag[X[k1+1]]=1; 
					}
					
					error1=constructionMisVal(DISS,TM1,X,W);
					if (error1==-1) {free(X);
						free(Flag);
						free(SumW);
						free(Sum);    
						for (i1=0;i1<=n;i1++)   
							{ free(TM1[i1]); free(TMopt[i1]); }    
					free(TM1); free(TMopt); return -1;}
					approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
					if (error1==-1) {free(X);
						free(Flag);
						free(SumW);
						free(Sum);    
						for (i1=0;i1<=n;i1++)   
							{ free(TM1[i1]); free(TMopt[i1]); }    
					free(TM1); free(TMopt); return -1;} 
					
					EQ=0.0;
					for (i1=1;i1<=n;i1++)
					{
						for (j1=i1+1;j1<=n;j1++)
							EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
					} 
					
					if ((a==0)||(EQ<EQmin))
					{
						a=1;
						EQmin=EQ;
						iopt=i;
						jopt=j;
						
						for (i1=1;i1<=n;i1++){
							for (j1=i1+1;j1<=n;j1++)
								TMopt[i1][j1]=TMopt[j1][i1]=TM[i1][j1];
						}			
					}         
					if (Option==2) break; 
				}   
			}
		}
	
		for (i1=1;i1<=n;i1++){
			for (j1=i1+1;j1<=n;j1++)
				TM[i1][j1]=TM[j1][i1]=TMopt[i1][j1];
			TM[i1][i1]=0.0;
		}
	}
	
	/* Cas de plusieurs processeurs */
	else{
		for (i=1;i<=n;i++){
			if(i%numprocs == myid ){
				for (j=i+1;j<=n;j++){
					if (W[i][j]!=0.0)
					{ 
						X[1]=i;
						X[2]=j;                    
						
						for (k=1;k<=n;k++)
							Flag[k]=0;
						Flag[i]=Flag[j]=1;
						for (k=1;k<=n;k++)
						{ 
							if (Flag[k]!=1)   
							{
								SumW[k]=W[k][i];
							}
						}
						
						for (k1=2;k1<=n-1;k1++)
						{
							MAX=MAX1=0.0;
							for (k=1;k<=n;k++)
							{ 
								if (Flag[k]!=1)   
								{             
									SumW[k]=SumW[k]+W[k][X[k1]];
									if ((MAX<SumW[k])||((MAX==SumW[k])&&(MAX1<Sum[k]))){ MAX=SumW[k]; X[k1+1]=k; MAX1=Sum[k];}    
								}              
							}
							Flag[X[k1+1]]=1; 
						}
						
						error1=constructionMisVal(DISS,TM1,X,W);
						if (error1==-1) {free(X);
							free(Flag);
							free(SumW);
							free(Sum);    
							for (i1=0;i1<=n;i1++)   
								{ free(TM1[i1]); free(TMopt[i1]); }    
						free(TM1); free(TMopt); return -1;}
						approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
						if (error1==-1) {free(X);
							free(Flag);
							free(SumW);
							free(Sum);    
							for (i1=0;i1<=n;i1++)   
								{ free(TM1[i1]); free(TMopt[i1]); }    
						free(TM1); free(TMopt); return -1;} 
						
						EQ=0.0;
						for (i1=1;i1<=n;i1++)
						{
							for (j1=i1+1;j1<=n;j1++)
								EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
						} 
						
						if ((a==0)||(EQ<EQmin))
						{
							a=1;EQmin=EQ;iopt=i;jopt=j;Teqmin[myid]=EQ;Tiopt[myid]=i;Tjopt[myid]=j;Ta[myid]=1;           
							
							for (i1=1;i1<=n;i1++)
							{
								for (j1=i1+1;j1<=n;j1++)
									TMopt[i1][j1]=TMopt[j1][i1]=TM[i1][j1];
							}			
						}         
						if (Option==2) break; 
					}   
				}
			}
		} 
		/* Syncrhonisation des resultats */
		if(myid>0){
			MPI_Send( &EQmin, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
			MPI_Send( &a, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
			MPI_Send( &iopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
			MPI_Send( &jopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );

			MPI_Recv( &leproc, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
			if(myid==leproc){
				for (i1=1;i1<=n;i1++){
					for (j1=i1+1;j1<=n;j1++)
						MPI_Send( &TMopt[i1][j1], 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
				}

			}
			for (i1=1;i1<=n;i1++){
				for (j1=i1+1;j1<=n;j1++){	
					MPI_Recv( &TM[i1][j1], 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &status );
					TM[j1][i1]=TM[i1][j1];
				}
				TM[i1][i1]=0.0;
			}
			if(myid<numprocs-1){
				for (i1=1;i1<=n;i1++){
					for (j1=i1+1;j1<=n;j1++){
						MPI_Send( &TM[i1][j1], 1, MPI_DOUBLE,myid+1, 0, MPI_COMM_WORLD );
					}
				}
			}
		}
		else{
			for(i=1; i<numprocs;i++){
				MPI_Recv( &Teqmin[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Ta[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tiopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tjopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
			}
			if(a==0)
				EQmin=LONG_MAX;
			leproc=0;
			for(i=1; i<numprocs;i++){
				if ( (Ta[i]==1) &&((Tiopt[i]==1)&&(Tjopt[i]==2))||(Teqmin[i]<=EQmin)){
					if( (   (Teqmin[i]==EQmin) && (Tiopt[i]<iopt)) || ((Teqmin[i]==EQmin) && (Tiopt[i]==iopt) && (Tjopt[i]<jopt)))  {
							EQmin=Teqmin[i];iopt=Tiopt[i];jopt=Tjopt[i];leproc=i;    
					}
					if( (Teqmin[i]<EQmin)){
							EQmin=Teqmin[i];iopt=Tiopt[i];jopt=Tjopt[i]; leproc=i;  
					}
				}
			}
			for(i=1; i<numprocs;i++)
				MPI_Send( &leproc, 1, MPI_INT,i, 0, MPI_COMM_WORLD );

			if(leproc!=myid){
				for (i1=1;i1<=n;i1++){
					for (j1=i1+1;j1<=n;j1++){	
						MPI_Recv(&TM[i1][j1], 1, MPI_DOUBLE, leproc, 0, MPI_COMM_WORLD, &status );
						TM[j1][i1]=TM[i1][j1];
					}
					TM[i1][i1]=0.0;
				}
				
			}
			for (i1=1;i1<=n;i1++){
				for (j1=i1+1;j1<=n;j1++){
					MPI_Send( &TM[i1][j1], 1, MPI_DOUBLE,myid+1, 0, MPI_COMM_WORLD );
				}
			}
		}
	}
	/*for (i1=1;i1<=n;i1++){
		for (j1=i1+1;j1<=n;j1++)
			printf ("\n%lf et mon myid = %d et iopt=%d et jopt=%d et EQmin=%d", TM[i1][j1],myid, iopt, jopt, EQmin);
	}*/

	free(X);
	free(Flag);
	free(SumW);
	free(Sum);    
	for (i=0;i<=n;i++)   
		{ free(TM1[i]); free(TMopt[i]); }  
	free(TM1);
	free(TMopt);
	
	printf("\n FIN");
	return 1;
}

/* MW-modified (part2) */
int parcour211(double **DISS,double **W,double **TM,int *Iternum,long int *ARETE, double *LONGUEUR)
{
	int i,j,i1,j1,iopt,jopt,*X,Iternumber,Option;
	double **TM1;
	double *Teqmin; /* Tableau des EQmin necessaire à la parallelisation */
	int *Tiopt, *Tjopt; /* Tableau des iopt et jopt necessaire à la parallelisation */
	double EQ,EQmin;
	
	int error1=1;
	TM1=(double **) malloc((n+1)*sizeof(double*));
	X=(int *) malloc((n+1)*sizeof(int)); 

	/* Allocation juste dans le cas d'une execution en parallele*/
	if(numprocs>1){
		Tiopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Tjopt=(int *) malloc((numprocs+1)*sizeof(int)); 
		Teqmin =(double*)malloc((numprocs+1)*sizeof(double)); 
	}

	EQmin=(double)LONG_MAX;


	for (i=0;i<=n;i++)
	{
		TM1[i]=(double*)malloc((n+1)*sizeof(double));    
		if (TM1[i]==NULL)
		{
			printf("Data matrix is too large"); 
			exit(5); 
		}  
	}    
	
	Iternumber=*Iternum;
	Option=floor1(DISS[0][0]);
	X[0]=1;
	
	Option=1;
	/* Cas d'un processeur */
	if(numprocs==1){
		for (i=1;i<=n;i++){  
			for (j=i+1;j<=n;j++)   {
				X[1]=i;
				X[2]=j;          
				construction(DISS,TM1,X,W); 
				if (error1==-1) {free(X);  
					for (i1=0;i1<=n;i1++)   
						free(TM1[i1]);   
				free(TM1); return -1;}            
				approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
				if (error1==-1) {
					free(X);  
					for (i1=0;i1<=n;i1++)   
						free(TM1[i1]);   
					free(TM1); return -1;}           
				
				EQ=0.0;
				for (i1=1;i1<=n;i1++)
				{
					for (j1=i1+1;j1<=n;j1++)
						EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
				} 
				
				if (((i==1)&&(j==2))||(EQ<EQmin))
				{
					EQmin=EQ;
					iopt=i;
					jopt=j;     
				}         
				if (Option==2) break;  
			}
		}
	}
	/* Cas de plusieurs processeurs */
	else{
		for (i=1;i<=n;i++){  
			if(i%numprocs == myid ){
				for (j=i+1;j<=n;j++){
					X[1]=i;
					X[2]=j;          
					construction(DISS,TM1,X,W); 
					if (error1==-1) {free(X);  
						for (i1=0;i1<=n;i1++)   
							free(TM1[i1]);   
					free(TM1); return -1;}            
					approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
					if (error1==-1) {free(X);  
						for (i1=0;i1<=n;i1++)   
							free(TM1[i1]);   
					free(TM1); return -1;}           
					
					EQ=0.0;
					for (i1=1;i1<=n;i1++)
					{
						for (j1=i1+1;j1<=n;j1++)
							EQ=EQ+W[i1][j1]*(DISS[i1][j1]-TM[i1][j1])*(DISS[i1][j1]-TM[i1][j1]);
					} 
					
					if (((i==1)&&(j==2))||(EQ<EQmin))
					{
						EQmin=EQ;iopt=i;jopt=j;Teqmin[myid]=EQ;Tiopt[myid]=i;Tjopt[myid]=j;            
					}         
					if (Option==2) break;  
				}
			}
		}	
		/* Synchronisation des resultats */       
		if(myid>0){
			MPI_Send( &EQmin, 1, MPI_DOUBLE,0, 0, MPI_COMM_WORLD );
			MPI_Send( &iopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );
			MPI_Send( &jopt, 1, MPI_INT,0, 0, MPI_COMM_WORLD );	
			MPI_Recv( &iopt, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
			MPI_Recv( &jopt, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status );
			if(myid<numprocs-1){
				MPI_Send( &iopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
				MPI_Send( &jopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			}
		}
		else{
			for(i=1; i<numprocs;i++){
				MPI_Recv( &Teqmin[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tiopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( &Tjopt[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );		
			}
			for(i=1; i<numprocs;i++){
				if (((Tiopt[i]==1)&&(Tjopt[i]==2))||(Teqmin[i]<EQmin)){
						EQmin=Teqmin[i];iopt=Tiopt[i];jopt=Tjopt[i];    
					}
			}
			MPI_Send( &iopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
			MPI_Send( &jopt, 1, MPI_INT,myid+1, 0, MPI_COMM_WORLD );
		}
	}



	if (n>10) Iternumber=2*n-3;
	else Iternumber=15;
	X[1]=iopt;
	X[2]=jopt;  
	
	construction(DISS,TM1,X,W);
	if (error1==-1) {free(X);  
		for (i=0;i<=n;i++)   
			free(TM1[i]);   
	free(TM1); return -1;} 
	approx_arb2 (DISS,TM1,TM,W,&Iternumber,ARETE,LONGUEUR);
	if (error1==-1) {free(X);  
		for (i=0;i<=n;i++)   
			free(TM1[i]);   
	free(TM1); return -1;}     
        
	free(X);  
	for (i=0;i<=n;i++)   
		free(TM1[i]);   
	free(TM1);
	return 1;
}


/* Tree reconstruction from incomplete matrix by MW-modified (part3) */
int constructionMisVal(double **D,double **DA,int *X,double **W)
{
	int i,j,i1,p,P,xi,a,*Y,NV,NR,PP,PN,option,*ARETE;
	double *A,am,c,am1,c1,k,k1,k2,k5,k6,k7,ki,**G1,S1,S2,
	S3,M,MR,M1,DIS,DIS1,Su,Sv,*L,**L1,**L2,
	**L3,*DIST,*DIST1;
	double Dummy1;
	
	
	Y=(int *) malloc((n+1)*sizeof(int));
	A=(double *) malloc((8*n+8)*sizeof(double));
	G1=(double **) malloc((2*n-2)*sizeof(double*));
	
	L=(double *) malloc((n+1)*sizeof(double));
	L1=(double **) malloc((2*n-2)*sizeof(double*));
	L2=(double **) malloc((2*n-2)*sizeof(double*));
	L3=(double **) malloc((2*n-2)*sizeof(double*));
	DIST=(double *) malloc(3*sizeof(double));
	DIST1=(double *) malloc(3*sizeof(double));
	ARETE=(int *) malloc(3*sizeof(int));
	
	for (i=0;i<=2*n-3;i++)
	{
		G1[i]=(double*)malloc((n+1)*sizeof(double));
		L1[i]=(double*)malloc((n+1)*sizeof(double));
		L2[i]=(double*)malloc((n+1)*sizeof(double));
		L3[i]=(double*)malloc((n+1)*sizeof(double)); 
		if ((G1[i]==NULL)||(L1[i]==NULL)||(L2[i]==NULL)||(L3[i]==NULL))
		{
			printf("Data matrix is too large"); 
			exit(5); 
		}
	}      
	
	Dummy1=100;
	if (D[X[1]][X[2]]<(n-1)*0.0002){ Dummy1=D[X[1]][X[2]]; D[X[1]][X[2]]=(n-1)*0.0002;}
	
	option=X[0];
	
	MR=0;
	A[1]=X[1];
	A[2]=X[2];
	A[3]=0;
	
	A[4]=D[X[1]][X[2]]; 
	
	for (i=1;i<=n;i++)
		Y[i]=1;
	Y[X[1]]=0;
	Y[X[2]]=0;
	P=1;
	DA[X[1]][X[1]]=0;
	DA[X[2]][X[2]]=0;
	DA[X[1]][X[2]]=A[4];
	DA[X[2]][X[1]]=A[4];
	
	xi=X[2];
	
	for (j=1;j<=n;j++)
	{
		if ((j!=X[1])&&(j!=X[2]))
		{
			k1=DA[X[1]][xi]-D[xi][j];
			k2=D[X[1]][j];
			L[j]=W[X[1]][j]+W[X[2]][j];
			L1[1][j]=2*(k1*W[X[2]][j]-k2*W[X[1]][j]);
			L2[1][j]=2*(-k1*W[X[2]][j]-k2*W[X[1]][j]);
			L3[1][j]=2*(W[X[1]][j]-W[X[2]][j]);
			G1[1][j]=k1*k1*W[X[2]][j]+k2*k2*W[X[1]][j];
		}
	}
	
	for (i=2;i<=n-1;i++)
	{    
		a=0;   
		if ((option!=2)||(i>3)) 
		{
			
			for (p=1;p<=P;p++)
			{
				for (j=1;j<=n;j++)
				{
					if ((Y[j]==1)&&(X[i+1]==j))  
					{
						xi=floor1(A[4*p-2]);
						k1=DA[X[1]][xi]-D[xi][j];
						k2=D[X[1]][j];
						Su=A[4*p-1];
						Sv=Su+A[4*p];
						
						k6=L1[p][j];
						k5=L2[p][j];
						k7=L3[p][j];
						k=L[j];
						M=G1[p][j]+MR+1.0;
						
						if ((2*k*Sv+k5<=0.00001)&&(k6+k7*Sv>=-0.00001))
						{ M=k*Sv*Sv+k5*Sv;
							am=Sv;
							c=0;
						}
						if (((c1=-(k6+k7*Sv)/(2*k))>=-0.00001)&&(-k5-2*k*Sv-k7*c1>=-0.00001)&&((M1=k*(Sv*Sv+c1*c1)+k7*Sv*c1+k5*Sv+c1*k6)<M))
						{ M=M1;
							am=Sv;
							c=c1;
						}
						if (((am1=-k5/(2*k))>=Su-0.00001)&&(-k6-k7*am1<=0.00001)&&((M1=k*am1*am1+k5*am1)<M)&&(am1<=Sv+0.00001))
						{ M=M1;
							c=0;
							am=am1;
						}
						if (((4*k*k-k7*k7)!=0)&&((am1=(-2*k*k5+k7*k6)/(4*k*k-k7*k7))<=Sv+0.00001)&&((c1=-(am1*k7+k6)/(2*k))>=-0.00001)&&((M1=k*(am1*am1+c1*c1)+k7*am1*c1+k5*am1+k6*c1)<M)&&(am1>=Su-0.00001))
						{
							M=M1;
							c=c1;
							am=am1;
						}
						
						if (((4*k*k-k7*k7)==0)&&(k5==k6)&&((-k5/(2*k))>=Su-0.00001)&&((c1=-k5/(2*k)-Su)>=-0.00001)&&((M1=k*(Su+c1)*(Su+c1)+k5*(Su+c1))<M))
						{
							M=M1;
							c=c1;
							am=Su;
						}
						
						if ((2*k==k7)&&(k5==k6)&&(k5>=-0.00001)&&((M1=k*Su*Su+k5*Su)<M))
						{
							M=M1;
							c=0;
							am=Su;
						}
						
						
						if ((2*k*Su+k5>=-0.00001)&&(-k6-k7*Su<=0.00001)&&((M1=k*Su*Su+k5*Su)<M))
						{ c=0;
							am=Su;
							M=M1;
						}
						
						if (((c1=-(Su*k7+k6)/(2*k))>=-0.00001)&&((k5+2*k*Su+c1*k7)>=-0.00001)&&((M1=k*Su*Su+k5*Su+c1*c1*k+c1*(Su*k7+k6))<M))
						{ 
							am=Su;
							M=M1;
							c=c1;
						}
						if ((a==0)||(MR>M+G1[p][j]))
						{
							a=1;
							MR=M+G1[p][j];
							DIS=am;
							DIS1=c;
							
							
							if (c<=(n-i+1.0)*0.00002) DIS1=(n-i+1.0)*0.00002;
							if (fabs(Sv-Su)<=2*(n-i+1)*0.00002) DIS=Sv-(Sv-Su)/(n-i+1.0);
							else 
							{   
								if (fabs(am-Su)<=(n-i+1.0)*0.00002) DIS=Su+(n-i+1.0)*0.00002;
								if (fabs(am-Sv)<=(n-i+1.0)*0.00002) DIS=Sv-(n-i+1.0)*0.00002;
							} 
							
							
							NV=j;
							NR=p;
						}
					}
				}
			}
			
		}
		
		if ((option==2)&&(i<=3))
		{
			NV=X[i+1];
			NR=ARETE[i-1];
			DIS=DIST[i-1];
			DIS1=DIST1[i-1];
		}  
		
		
		/*            Updating data in the tables          */
		p=NR;
		PP=P;
		if (fabs(DIS-A[4*p-1])<=0.00001)
			DIS=A[4*p-1];
		if (fabs(DIS-A[4*p-1]-A[4*p])<=0.00001)
		{
			DIS=A[4*p-1]+A[4*p];
		}
		X[i+1]=NV;
		xi=floor1(A[4*p-2]);
		Y[NV]=0;
		if ((DIS-A[4*p-1]>0)&&(DIS-A[4*p-1]-A[4*p]<0))
		{
			c=A[4*p];
			A[4*p]=DIS-A[4*p-1];
			P=P+1;
			A[4*P-3]=P;
			A[4*P-2]=A[4*p-2];
			A[4*P-1]=DIS;
			A[4*P]=A[4*p-1]+c-DIS;
		}
		if (DIS1>0)
		{
			P=P+1;
			A[4*P-3]=P;
			A[4*P-2]=X[i+1];
			A[4*P-1]=DIS;
			A[4*P]=DIS1;
		}
		DA[X[1]][X[i+1]]=DIS+DIS1;
		DA[X[i+1]][X[1]]=DIS+DIS1;
		
		
		/* The induction formula */
		
		for (j=2;j<=i;j++)
		{
			
			if (((DA[X[1]][X[j]]+DA[X[1]][xi]-DA[xi][X[j]])/2)<=DIS)
				DA[X[j]][X[i+1]]=DIS+DIS1+DA[X[j]][xi]-DA[X[1]][xi];
			
			else
				DA[X[j]][X[i+1]]=DIS1-DIS+DA[X[1]][X[j]];
			
			
			DA[X[i+1]][X[j]]=DA[X[j]][X[i+1]];
			DA[X[i+1]][X[i+1]]=0;
			
		}
		
		if ((P==PP+2)||((P==PP+1)&&(DIS1==0)))
			PN=PP+1;
		else
			PN=PP;
		if ((DIS-A[4*p-1]>0)&&(DIS-A[4*p-1]-c<0))
		{
			for(j=1;j<=n;j++)
			{
				if (Y[j]==1)
				{
					L1[PP+1][j]=L1[p][j];
					L2[PP+1][j]=L2[p][j];
					L3[PP+1][j]=L3[p][j];
					G1[PP+1][j]=G1[p][j];
				}
			}
		}
		
		/* BLOCK 1 */
		for (p=1;p<=PN;p++)
		{
			xi=floor1(A[4*p-2]);
			for (j=1;j<=n;j++)
			{
				if (Y[j]==1)
				{ 
					if (((DA[X[1]][xi]+DA[X[1]][X[i+1]]-DA[X[i+1]][xi])/2)<=A[4*p-1]+0.00001)
					{
						ki=D[X[i+1]][j]+DA[X[1]][xi]-DA[X[i+1]][xi];
						L1[p][j]=L1[p][j]-2*ki*W[X[i+1]][j];
						L2[p][j]=L2[p][j]-2*ki*W[X[i+1]][j];
						L3[p][j]=L3[p][j]+2*W[X[i+1]][j];
					}
					else
					{
						ki=D[X[i+1]][j]-DA[X[i+1]][X[1]];
						L1[p][j]=L1[p][j]-2*ki*W[X[i+1]][j];
						L2[p][j]=L2[p][j]+2*ki*W[X[i+1]][j];
						L3[p][j]=L3[p][j]-2*W[X[i+1]][j];
					}
					
					if (p==1)
						L[j]=L[j]+W[X[i+1]][j];
					G1[p][j]=G1[p][j]+ki*ki*W[X[i+1]][j];
				}
			}
		}
		
		/* BLOCK 2 */
		if (DIS1>0)
		{
			for (j=1;j<=n;j++)
			{
				if (Y[j]==1)
				{
					k1=DA[X[1]][X[i+1]]-D[X[i+1]][j];
					k2=D[X[1]][j];
					S1=W[X[1]][j]+W[X[i+1]][j];
					S2=k2*W[X[1]][j];
					S3=k1*k1*W[X[i+1]][j]+k2*k2*W[X[1]][j];
					for (i1=2;i1<=i;i1++)
					{
						ki=D[X[i1]][j]+DA[X[1]][X[i+1]]-DA[X[i+1]][X[i1]];
						S1=S1+W[X[i1]][j];
						S2=S2+ki*W[X[i1]][j];
						S3=S3+ki*ki*W[X[i1]][j];
					}
					L[j]=S1;
					L1[P][j]=2*(k1*W[X[i+1]][j]-S2);
					L2[P][j]=2*(-k1*W[X[i+1]][j]-S2);
					L3[P][j]=2*(S1-2*W[X[i+1]][j]);
					G1[P][j]=S3;
				}
			}
		}
	}
	
	
	if (Dummy1<(n-1)*0.0002) D[X[1]][X[2]]=Dummy1;
	
	free(L);
	free(Y);
	free(A);
	free(DIST);
	free(DIST1);
	free(ARETE);
	
	for (i=0;i<=2*n-3;i++)
	{ 
		free(G1[i]);
		free(L1[i]);
		free(L2[i]);
		free(L3[i]);
	}    
	free(G1);
	free(L1);
	free(L2);
	free(L3);
	
	return 1;
}


int tryout(double **DI, double **W)
{
	
	int i,j,k,l;
	double S, S1, MAX, MIN=1000.0, **DIN, **WN;
	
	int Flag,Flag1,Flag1old=0,tour=0;
	char method ='A';
	double MissPerc;
	
	DIN = (double **)malloc((n+1)*sizeof(double*)); 
	WN = (double **)malloc((n+1)*sizeof(double*)); 
	
	for (i=0; i<=n; i++)
	{	DIN[i]=(double *)malloc((n+1)*sizeof(double)); 
		WN[i]=(double *)malloc((n+1)*sizeof(double));
	}
	
	Label: 
	tour++; Flag1=0;
	for (i=1; i<=n; i++)
	{
		for (j=1; j<=n; j++)
		{  		  		
			DIN[i][j]=DI[j][i];
			WN[i][j]=W[j][i];
			if (DI[i][j]==-99.00) { Flag1++; }
		}
	}
	
	if (tour==1) MissPerc = (n*(n-1.0)-Flag1)/(n*(n-1.0));
	
	for (i=1; i<=n; i++)
	{
		for (j=i+1; j<=n; j++)
		{  		  		
			if ((DI[i][j]==-99.00)&&(DI[j][i]==-99.00)) 
			{
				S=0.0; S1=0.0;
				
				if (n > (3.0 + 2.0/(MissPerc*MissPerc*MissPerc))) /*Additive procedure*/
				{
					method ='A';
					for (k=1; k<=n; k++)
					{
						if ((k!=i)&&(k!=j))
						{
							for (l=k+1; l<=n; l++)
							{  					
								if ((l!=i)&&(l!=j))
								{
									if ((W[k][l]!=0.0)&&(W[i][k]!=0.0)&&(W[j][l]!=0.0)&&(W[i][l]!=0.0)&&(W[j][k]!=0.0)&&(fabs(DI[i][k]+DI[j][l]-DI[i][l]-DI[j][k]) > 0.00001))
									{
										MAX = (DI[i][k]+DI[j][l])>(DI[i][l]+DI[j][k]) ? (DI[i][k]+DI[j][l]) : (DI[i][l]+DI[j][k]);												
										if ((S==0)||(MIN > MAX - DI[k][l])) MIN = S = MAX - DI[k][l];
										S1 = S1 + 1;
									}
								}
							}
						}
					}
				}
				else /*Ultrametric procedure */
				{ 
					method ='U'; S=0.0; S1=0.0;
					for (k=1; k<=n; k++)
					{
						if ((k!=i)&&(k!=j))
						{
							if ((W[i][k]!=0.0)&&(W[j][k]!=0.0)&&(fabs(DI[i][k]-DI[j][k]) > 0.00001))
							{
								MAX = DI[i][k] > DI[j][k] ? DI[i][k] : DI[j][k];												
								if ((S==0)||(MIN > MAX)) MIN = S = MAX;
								S1 = S1 + 1;
							}
							
						}
					}		
				}
				
				if (S<=0.0) S = 0.0005;
				if ((S1!=0.0)&&(method =='A')) { DIN[i][j]=DIN[j][i]=S; WN[i][j]=WN[j][i]=0.5; }
				if ((S1!=0.0)&&(method =='U')) { DIN[i][j]=DIN[j][i]=S; WN[i][j]=WN[j][i]=0.5; }
				
			}
		}   
	}
	
	Flag1=Flag=0;
	for (i=1; i<=n; i++)
	{
		for (j=1; j<=n; j++)
		{  		  				
			DI[i][j]=DIN[i][j];
			W[i][j]=WN[i][j];
			if (DI[i][j]==-99.00) { Flag=1; Flag1++; }
		}
	}
	
	if ((Flag==1)&&(Flag1!=Flag1old)) 
		{ if (Flag1old!=0) Flag1old=Flag1; goto Label; }
	
	for (i=0; i<=n; i++)
		{	free(DIN[i]); free(WN[i]); }
	
	free(DIN);
	free(WN);
	
	return Flag;
	
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         BIONJ program                                     ;
;                                                                           ;
;                         Olivier Gascuel                                   ;
;                                                                           ;                                                                        ;
;                         LIRMM - Montpellier- France                       ;
;                         gascuel@lirmm.fr                                  ;
;                                                                                      ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

/*;;;;;;;;;;;  INPUT, OUTPUT, INITIALIZATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                                                                           ;
;              The delta matrix is read from the input-file.                ;
;              It is recommended to put it and the executable in            ;
;              a special directory. The input-file and output-file          ;
;              can be given as arguments to the executable by               ;
;              typing them after the executable (Bionj input-file           ;
;              output-file) or by typing them when asked by the             ;
;              program. The input-file has to be formated according         ;
;              the PHYLIP standard. The output file is formated             ;
;              according to the NEWSWICK standard.                          ;
;                                                                           ;
;              The lower-half of the delta matrix is occupied by            ;
;              dissimilarities. The upper-half of the matrix is             ;
;              occupied by variances. The first column                      ;
;              is initialized as 0; during the algorithm some               ;
;              indices are no more used, and the corresponding              ;
;              positions in the first column are set to 1.                  ;
;                                                                           ;
;              This delta matix is made symmetrical using the rule:         ;
;              Dij = Dji <- (Dij + Dji)/2. The diagonal is set to 0;        ;
;              during the further steps of the algorithm, it is used        ;
;              to store the sums Sx.                                        ;
;                                                                           ;
;              A second array, trees, is used to store taxon names.         ;
;              During the further steps of the algoritm, some               ;
;              positions in this array are emptied while the others         ;
;              are used to store subtrees.                                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


/*;;;;;;;;;;;;;;;;;;;;;;;;;; Initialize        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function reads an input file and return the            ;
;               delta matrix and trees: the list of taxa.                   ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              FILE *input    : pointer to input file                       ;
;              int n          : number of taxa                              ;
;              char **trees   : list of taxa                                ;
;                                                                           ;
; return value:                                                             ;
;              float **delta : delta matrix                                 ;
;              char *trees    : list of taxa                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Initialize(float **delta, FILE *input, int n, POINTERS *trees)
{
	int lig;                                          /* matrix line       */
	int col;                                          /* matrix column     */
	float distance;
	char name_taxon[LEN];                             /* taxon’s name      */
	WORD1 *name;
	
	for(lig=1; lig <= n; lig++)
	{
		fscanf(input,"%s",name_taxon);                  /* read taxon’s name */
		name=(WORD1 *)calloc(1,sizeof(WORD1));            /* taxon’s name is   */
		if(name == NULL)                                /* put in trees      */
		{
			printf("Out of memories !!");
			exit(0);
		}
		else
		{
			strcpy(name->name,name_taxon);
			name->suiv=NULL;
			trees[lig].head=name;
			trees[lig].tail=name;
			for(col= 1; col <= n; col++)
			{
				fscanf(input,"%f",&distance);             /* read the distance  */
				delta[lig][col]=distance;
			}
		}
	}
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Print_output;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function prints out the tree in the output file.       ;
;                                                                           ;
; input       :                                                             ;
;              POINTERS *trees : pointer to the subtrees.                   ;
;              int i          : indicate the subtree i to be printed.       ;
:              FILE *output   : pointer to the output file.                 ;
;                                                                           ;
; return value: The phylogenetic tree in the output file.                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


void Print_output(int i, POINTERS *trees, FILE *output)
{
	WORD1 *parcour;
	parcour=trees[i].head;
	while(parcour != NULL)
	{
		fprintf(output,"%s",parcour->name);
		parcour=parcour->suiv;
	}
	
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         Main program                                      ;
;                                                                           ;
;                         argc is the number of arguments                   ;
;                         **argv contains the arguments:                    ;
;                         the first argument has to be BIONJ;               ;
;                         the second is the inptu-file;                     ;
;                         the third is the output-file.                     ;
;                         When the input and output files are               ;
;                         not given, the user is asked for them.            ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


void BioNJ_main (double **DI, double **DA)
{
	FILE *input=NULL;                            /* pointer to input file       */
	FILE *output=NULL;                           /* pointer to output file      */
	POINTERS *trees=NULL;                        /* list of subtrees            */
	char *Name_fich1;                       /* name of the input file      */
	char *Name_fich2;                       /* name of the output file     */
	char *chain1;                           /* stringized branch-length    */
	char *chain2;                           /* idem                        */
	int *a, *b;                             /* pair to be agglomerated     */
	float **delta;                          /* delta matrix                */
	float la;                               /* first taxon’s branch-length */
	float lb;                               /* second taxon’s branch-length*/
	float vab;                              /* variance of Dab             */
	float lamda;
	int i;
	int ok;
	int r;                                  /* number of subtrees          */
	//int n;                                  /* number of taxa              */
	int x, y;
	clock_t clock_start, clock_end;
	double t;
	int Pr,j;
	int last1, last2, last3;
	
	double *T1,L,Lii,Ljj,l1,l2,l3;
	int *T,ii,jj,n1;
	
	T1=(double *) malloc((n+1)*sizeof(double));
	T=(int *) malloc((n+1)*sizeof(int));

	for (i=1;i<=n;i++)
	{
		T[i]=i;
		T1[i]=0.0;
	}
	
	/*   Allocation of memories    */
	
	a=(int*)calloc(1,sizeof(int));
	b=(int*)calloc(1,sizeof(int));
	
	/*      Create the delta matrix     */
	
	delta=(float **)calloc(n+1,sizeof(float*));
	for(i=1; i<= n; i++)
	{
		delta[i]=(float *)calloc(n+1, sizeof(float));
		if(delta[i] == NULL)
		{
			{ printf("Data matrix is too large - BioNJ"); /*m_error=-1;*/ return;}
		}
	}
	
	for(i=1; i <= n; i++)
	{
		for(j=1; j <= n; j++)
			delta[i][j]=DI[i][j];     
	}
	
	for(i=1; i<=n; i++)       	          			                
		delta[i][0]=0.0;		          
	
	
	/*   initialise and symmetrize the running delta matrix    */
	r=n;
	*a=0;
	*b=0;
	
	while (r > 3)                             /* until r=3                 */
	{
		Compute_sums_Sx(delta, n);             /* compute the sum Sx       */
		Best_pair(delta, r, a, b, n);          /* find the best pair by    */
		vab=Variance(*a, *b, delta);           /* minimizing (1)           */
		la=Branch_length(*a, *b, delta, r);    /* compute branch-lengths   */
		lb=Branch_length(*b, *a, delta, r);    /* using formula (2)        */
		lamda=Lamda(*a, *b, vab, delta, n, r); /* compute lambda* using (9)*/
		
		for(i=1; i <= n; i++)
		{
			if(!Emptied(i,delta) && (i != *a) && (i != *b))
			{
				if(*a > i)
				{
					x=*a;
					y=i;
				}
				else
				{
					x=i;
					y=*a;                           /* apply reduction formulae */
				}                                  /* 4 and 10 to delta        */
				delta[x][y]=Reduction4(*a, la, *b, lb, i, lamda, delta);
				delta[y][x]=Reduction10(*a, *b, i, lamda, vab, delta);
			}
		}
		
		/* Mise a jour de DA */
		if (la<0.00001) la=0.00001;
		if (lb<0.00001) lb=0.00001;
		n1=r; ii=*b; jj=*a; Lii=lb; Ljj=la;	  
		
		for (i=1;i<=n;i++)
		{ 
			if (T[i]==ii) T1[i]=T1[i]+Lii;
			if (T[i]==jj) T1[i]=T1[i]+Ljj;
		}
		for (j=1;j<=n;j++)
		{ 
			if (T[j]==jj)
			{
				for (i=1;i<=n;i++)
				{
					if (T[i]==ii)
					{
						DA[i][j]=T1[i]+T1[j];
						DA[j][i]=DA[i][j];
					}
				}
			}
		}
		
		for (j=1;j<=n;j++)
			if (T[j]==ii) T[j]=jj;
		/* end of bloc update DA */
		delta[*b][0]=1.0;                     /* make the b line empty     */
		r=r-1;                                /* decrease r                */
	}
	
	
	Finish(delta, n, trees, output,&l1,&l2,&l3,&last1,&last2,&last3);       /* compute the branch-lengths*/

	if (l1<0.00001) l1=0.00001;
	if (l2<0.00001) l2=0.00001;
	if (l3<0.00001) l3=0.00001;
	
	/*Il reste 3 sommets */
	for (j=1;j<=n;j++)
	{
		for (i=1;i<=n;i++)
		{
			if ((T[j]==last1)&&(T[i]==last2))
			{
				DA[i][j]=T1[i]+T1[j]+l1+l2;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==last1)&&(T[i]==last3))
			{
				DA[i][j]=T1[i]+T1[j]+l1+l3;
				DA[j][i]=DA[i][j];
			}
			if ((T[j]==last2)&&(T[i]==last3))
			{
				DA[i][j]=T1[i]+T1[j]+l2+l3;
				DA[j][i]=DA[i][j];
			}
		}
		DA[j][j]=0;
	}
	
	free(T);
	free(T1);
	/* end of bloc Il reste 3 sommets */
	
	free(delta);
	return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                             Utilities                                     ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/



/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Symmetrize  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifies if the delta matrix is symmetric;    ;
;               if not the matrix is made symmetric.                        ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              int n          : number of taxa                              ;
;                                                                           ;
; return value:                                                             ;
;              int symmetric  : indicate if the matrix has been made        ;
;                               symmetric or not                            ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Symmetrize(float **delta, int n)
{
	int lig;                                         /* matrix line        */
	int col;                                         /* matrix column      */
	float value;                                     /* symmetrized value  */
	int symmetric;
	
	symmetric=1;
	for(lig=1; lig  <=  n; lig++)
	{
		for(col=1; col< lig; col++)
		{
			if(delta[lig][col] != delta[col][lig])
			{
				value= (delta[lig][col]+delta[col][lig])/2;
				delta[lig][col]=value;
				delta[col][lig]=value;
				symmetric=0;
			}
		}
	}
	if(!symmetric)
		printf("The matrix is not symmetric");
	return(symmetric);
}




/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Concatenate ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function concatenates a string to another.             ;
;                                                                           ;
; input       :                                                             ;
;      char *chain1    : the string to be concatenated.                     ;
;      int ind         : indicate the subtree to which concatenate the      ;
;                        string                                             ;
;      POINTERS *trees  : pointer to subtrees.                              ;
;      int post        : position to which concatenate (front (0) or        ;
;                        end (1))                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post)
{
	WORD1 *bran;
	
	bran=(WORD1 *)calloc(1,sizeof(WORD1));
	if(bran == NULL)
	{
		printf("Out of memories");
		exit(0);
	}
	else
	{
		strcpy(bran->name,chain1);
		bran->suiv=NULL;
	}
	if(post == 0)
	{
		bran->suiv=trees[ind].head;
		trees[ind].head=bran;
	}
	else
	{
		trees[ind].tail->suiv=bran;
		trees[ind].tail=trees[ind].tail->suiv;
	}
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Distance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve ant return de distance between taxa  ;
;               i and j from the delta matrix.                              ;
;                                                                           ;
; input       :                                                             ;
;               int i          : taxon i                                    ;
;               int j          : taxon j                                    ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               float distance : dissimilarity between the two taxa         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Distance(int i, int j, float **delta)
{
	if(i > j)
		return(delta[i][j]);
	else
		return(delta[j][i]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Variance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve and return the variance of the       ;
;               distance between i and j, from the delta matrix.            ;
;                                                                           ;
; input       :                                                             ;
;               int i           : taxon i                                   ;
;               int j           : taxon j                                   ;
;               float **delta  : the delta matrix                           ;
;                                                                           ;
; return value:                                                             ;
;               float distance : the variance of  Dij                       ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Variance(int i, int j, float **delta)
{
	if(i > j)
		return(delta[j][i]);
	else
		return(delta[i][j]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Emptied ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifie if a line is emptied or not.          ;
;                                                                           ;
; input       :                                                             ;
;               int i          : subtree (or line) i                        ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               0              : if not emptied.                            ;
;               1              : if emptied.                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Emptied(int i, float **delta)      /* test if the ith line is emptied */
{
	return((int)delta[i][0]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Sum_S;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function retrieves the sum Sx from the diagonal       ;
;                of the delta matrix.                                       ;
;                                                                           ;
;  input       :                                                            ;
;               int i          : subtree i                                  ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
;  return value:                                                            ;
;                float delta[i][i] : sum Si                                 ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Sum_S(int i, float **delta)          /* get sum Si form the diagonal */
{
	return(delta[i][i]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;Compute_sums_Sx;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function computes the sums Sx and store them in the    ;
;               diagonal the delta matrix.                                  ;
;                                                                           ;
; input       :                                                             ;
;     	         float **delta : the delta matrix.                      ;
;     	         int n          : the number of taxa                    ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Compute_sums_Sx(float **delta, int n)
{
	float sum;
	int i;
	int j;
	
	for(i= 1; i <= n ; i++)
	{
		if(!Emptied(i,delta))
		{
			sum=0;
			for(j=1; j <=n; j++)
			{
				if(i != j && !Emptied(j,delta))           /* compute the sum Si */
					sum=sum + Distance(i,j,delta);
			}
		}
		delta[i][i]=sum;                           /* store the sum Si in */
	}                                               /* delta’s diagonal    */
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Best_pair;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function finds the best pair to be agglomerated by    ;
;                minimizing the agglomerative criterion (1).                ;
;                                                                           ;
;  input       :                                                            ;
;                float **delta : the delta matrix                           ;
;                int r          : number of subtrees                        ;
;                int *a         : contain the first taxon of the pair       ;
;                int *b         : contain the second taxon of the pair      ;
;                int n          : number of taxa                            ;
;                                                                           ;
;  return value:                                                            ;
;                int *a         : the first taxon of the pair               ;
;                int *b         : the second taxon of the pair              ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Best_pair(float **delta, int r, int *a, int *b, int n)
{
	float Qxy;                         /* value of the criterion calculated*/
	int x,y;                           /* the pair which is tested         */
	float Qmin;                        /* current minimun of the criterion */
	
	Qmin=1.0e300;
	for(x=1; x <= n; x++)
	{
		if(!Emptied(x,delta))
		{
			for(y=1; y < x; y++)
			{
				if(!Emptied(y,delta))
				{
					Qxy=Agglomerative_criterion(x,y,delta,r);
					if(Qxy < Qmin-0.000001)
					{
						Qmin=Qxy;
						*a=x;          
						*b=y;
					}
				}  
			}
		}
	}
}


/*;;;;;;;;;;;;;;;;;;;;;;Finish_branch_length;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description :  Compute the length of the branch attached                 ;
;                 to the subtree i, during the final step                   ;
;                                                                           ;
;  input       :                                                            ;
;                int i          : position of subtree i                     ;
;                int j          : position of subtree j                     ;
;                int k          : position of subtree k                     ;
;                float **delta :                                            ;
;                                                                           ;
;  return value:                                                            ;
;                float length  : The length of the branch                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Finish_branch_length(int i, int j, int k, float **delta)
{
	float length;
	length=0.5*(Distance(i,j,delta) + Distance(i,k,delta)
	-Distance(j,k,delta));
	return(length);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Finish;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function compute the length of the lasts three        ;
;                subtrees and write the tree in the output file.            ;
;                                                                           ;
;  input       :                                                            ;
;                float **delta  : the delta matrix                          ;
;                int n           : the number of taxa                       ;
;                WORD1 *trees   : list of subtrees                           ;
;                                                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Finish(float **delta, int n, POINTERS *trees, FILE *output, double *l1, double *l2, double *l3, int *last1, int *last2, int *last3)
{
	int l=1;
	int i=0;
	float length;
	char *str;
	WORD1 *bidon;
	WORD1 *ele;
	int last[3];                            /* the last three subtrees     */
	
	while(l <= n)
	{                                       /* find the last tree subtree  */
		if(!Emptied(l, delta))
		{
			last[i]=l;
			i++;
		}
		l++;
	}
	
	*last1 = last[0];
	*last2 = last[1];
	*last3 = last[2];
	
	*l1 = length=Finish_branch_length(last[0],last[1],last[2],delta);

	*l2 = length=Finish_branch_length(last[1],last[0],last[2],delta);
	
	*l3 = length=Finish_branch_length(last[2],last[1],last[0],delta);

	
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;   
;                          Formulae                                         ;
;                                                                           ;  
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


float Agglomerative_criterion(int i, int j, float **delta, int r)
{
	float Qij;
	Qij=(r-2)*Distance(i,j,delta)                           /* Formula (1) */
	-Sum_S(i,delta)
	-Sum_S(j,delta); 
	
	return(Qij);                       
}


float Branch_length(int a, int b, float **delta, int r)
{
	float length;
	length=0.5*(Distance(a,b,delta)                         /* Formula (2) */
	+(Sum_S(a,delta)
	-Sum_S(b,delta))/(r-2)); 
	return(length);                                   
}


float Reduction4(int a, float la, int b, float lb, int i, float lamda,float **delta)
{
	float Dui;
	Dui=lamda*(Distance(a,i,delta)-la)+(1-lamda)*(Distance(b,i,delta)-lb);                /* Formula (4) */
	return(Dui);
}


float Reduction10(int a, int b, int i, float lamda, float vab,float **delta)
{
	float Vci;
	Vci=lamda*Variance(a,i,delta)+(1-lamda)*Variance(b,i,delta)-lamda*(1-lamda)*vab;                              /*Formula (10)  */
	return(Vci);
}

float Lamda(int a, int b, float vab, float **delta, int n, int r)
{
	float lamda=0.0;
	int i;
	
	if(vab==0.0)
		lamda=0.5;
	else
	{
		for(i=1; i <= n ; i++)
		{
			if(a != i && b != i && !Emptied(i,delta))
				lamda=lamda + (Variance(b,i,delta) - Variance(a,i,delta));
		}
		lamda=0.5 + lamda/(2*(r-2)*vab);             
	}                                           /* Formula (9) and the  */
	if(lamda > 1.0)                               /* constraint that lamda*/
		lamda = 1.0;								/* belongs to [0,1]     */
	if(lamda < 0.0)
		lamda=0.0;
	return(lamda);
}

int ExtraireDonnees(const char * chaine, char *champs, char * contenu){
	
	int cpt=0,i;
	int egale=0;
	
	if(chaine[0] != '-')
		return 0;
	
	for(i=1;i<strlen(chaine);i++){
		
		if (chaine[i] == '='){
			egale = 1;
			champs[cpt] = '\0';
			cpt=0;
			continue;
		}
		if(egale)
			contenu[cpt++] = chaine[i];
		else
			champs[cpt++] = chaine[i];
	}
	contenu[cpt] = '\0';
	
	if (!egale)
		return 0;
	
	return 1;
}


/* Lecture des options passé par la ligne de commande */

void LireOptions(char **argv, int argc){
	
	int i;
	char * champs  = (char*)malloc(100);
	char * contenu = (char*)malloc(100);
	
	//=valeur par defaut
	p_optionTR = 1;
	p_optionMiss = 2;
	p_method = 2;
	p_optionFunction = 1;
	p_iternumber = 5;
	p_option = 1;
	p_option1 = 1;
	p_p = 0;
	p_k = 5;
	p_boot = 2;
	p_nbRep = 10;
	
	m_sequenceMethod = 32887;	
	m_gapPenalty=0;  
	m_PEMV = 0;	
	m_paramA = -1;
	m_dataType = 1;
	
	for(i=1;i<argc;i++){
		if(ExtraireDonnees(argv[i],champs,contenu)){
			
		
			//== method sequence ======================
			if(strcmp("sequenceMethod",champs) == 0){
				m_sequenceMethod = atoi(contenu);
			}
			//== gapPenalty ======================
			else if(strcmp("gapPenalty",champs) == 0){
				m_gapPenalty = atoi(contenu);
			}
			//== PEMV ======================
			else if(strcmp("PEMV",champs) == 0){
				m_PEMV = atoi(contenu);
			}
			//== k ======================
			else if(strcmp("paramA",champs) == 0){
				m_paramA = atoi(contenu);
			}
			//== nombre de replicats ======================
			else if(strcmp("nbRep",champs) == 0){
				p_nbRep = atoi(contenu);
			}
			//== type des données 1=distance 2=sequences ==
			else if(strcmp("dataType",champs) == 0){
				m_dataType = atoi(contenu);
			}
			//== output newick ==================== 
			else if(strcmp("newickFile",champs) == 0){
				strcpy(p_newick,contenu);
			}
			//== output bootstrap ==================== 
			else if(strcmp("outbootstrap",champs) == 0){
				strcpy(p_outboot,contenu);
			}
			//== output matrix ======================= 
			else if(strcmp("outmatrix",champs) == 0){
				strcpy(p_matrix,contenu);
			}
			//== input file ====================== 
			else if(strcmp("input",champs) == 0){
				strcpy(p_input,contenu);
			}
			//== output file ======================
			else if(strcmp("output",champs) == 0){
				strcpy(p_output,contenu);
			}
			//== statistique ======================
			else if(strcmp("stat",champs) == 0){
				strcpy(p_stat,contenu);
			}
			//== modele ======================
			else if(strcmp("modele",champs) == 0){
				strcpy(p_modele,contenu);
			}
			//== k ======================
			else if(strcmp("k",champs) == 0){
				p_k = atoi(contenu);
			}
			//== p ======================
			else if(strcmp("p",champs) == 0){
				p_p = atoi(contenu);
			}
			//== option1 ======================
			else if(strcmp("option1",champs) == 0){
				p_option1 = atoi(contenu);
			}
			//== option ======================
			else if(strcmp("option",champs) == 0){
				p_option = atoi(contenu);
			}
			//== iternumber ======================
			else if(strcmp("iternumber",champs) == 0){
				p_iternumber = atoi(contenu);
			}
			//== optionFunction ======================
			else if(strcmp("optionFunction",champs) == 0){
				p_optionFunction = atoi(contenu);
			}
			//== optionTR ======================
			else if(strcmp("optionTR",champs) == 0){
				p_optionTR = atoi(contenu);
			}
			//== optionMiss ======================
			else if(strcmp("optionMiss",champs) == 0){
				p_optionMiss = atoi(contenu);
			}
			//== method ======================
			else if(strcmp("method",champs) == 0){
				p_method = atoi(contenu);
			}
			//== bootstrap ======================
			else if(strcmp("bootstrap",champs) == 0){
				p_boot = atoi(contenu);
			}
			else{
				printf("Unknown parameter --%s--!!",champs);
				exit(8);
			}
		}		
	}
	free(champs);
	free(contenu);
}


int f_CalculerDistances(FILE *pfile, char *outFile,char **pmatNomsEspeces, double **pmatDistances, 
char **pmatSites, int pintNbreEspeces, int pintNbreSites, int pintIdMethod, 
double pdblPenalty, int pintPEMV, double pdblParamA,int ecrireResultats,int pblnWarningDejaAffiche)
{
	int intRetour=0;			/* Retour de fonction*/
	int blnNucleo = 1;	/* Type de données du fichier source (nucléotides/protéines)*/
	int intMethode;			/* Méthode de calcul*/
	
	intRetour = f_LireFichierSource(pfile, pmatNomsEspeces, pmatDistances,
	pmatSites, pintNbreEspeces, pintNbreSites,
	&blnNucleo, pintIdMethod,
	pblnWarningDejaAffiche);
	if (intRetour != 0)
	{	p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, 
		pintNbreEspeces, pintNbreSites);
		return(-1);
	}
	
	intMethode = pintIdMethod;
	switch (intMethode)
	{	
		case ID_SM_UNCOR_DIST : /* Uncorrected distances*/
							if (blnNucleo)
								p_UncorrectedDistancesNucleo(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites, pdblPenalty);
							else
								p_UncorrectedDistancesProt(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites, pdblPenalty);
							break;
		case ID_SM_JUKES_CANTOR : /* Jukes-Cantor*/
							if (blnNucleo)
							{	if (pintPEMV == 0)
								intRetour = Compute_VB(pmatSites, pmatDistances, pintNbreEspeces, 
								pintNbreSites, pdblPenalty);	
								else
									p_JukesCantorNucleo(pmatSites, pmatDistances, pintNbreEspeces, 
								pintNbreSites, pdblPenalty);
							}
							else
								p_JukesCantorProt(pmatSites, pmatDistances, pintNbreEspeces, 
							pintNbreSites, pdblPenalty);
							break;
		case ID_SM_TAJIMA_NEI : /* Tajima-Nei*/
							if (blnNucleo)
								p_TajimaNei(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites);
							else
							{	printf("This method is not applicable with proteins");
								p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
							return(-1);	}
							break;
		case ID_SM_KIMURA_2_PARAM : /* Kimura 2-Parameters*/
							if (blnNucleo)
							{	if (pintPEMV == 0)
								intRetour = Compute_VB_K2P(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites);
								else
									p_Kimura2Parameter(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites);
							}
							else
							{	printf("This method is not applicable with proteins");
								p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
							return(-1);	}
							break;
		case ID_SM_TAMURA : /* Tamura*/
							if (blnNucleo)
								p_Tamura(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites);
							else
							{	printf("This method is not applicable with proteins");
								p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
							return(-1);	}
							break;
		case ID_SM_JIN_NEI_GAMMA : /* Jin-Nei Gamma*/
							if (blnNucleo)
							{	intRetour = f_JinNeiGamma(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites, 
								pdblParamA);
								if (intRetour != 0)
								{	p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
									return(-1);
								}
							}
							else
							{	printf("This method is not applicable with proteins");
								p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
							return(-1);	}
							break;
		case ID_SM_LOGDET : /* LogDet*/
							if (blnNucleo)
								p_LogDet(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites);
							else
							{	printf("This method is not applicable with proteins");
								p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
							return(-1);	}
							break;
		case ID_SM_F84 : /* F84 */
							if (blnNucleo)
								p_F84(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites);
							else
							{	printf("This method is not applicable with proteins");
								p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
							return(-1);	}
							break;
		case ID_SM_KIMURA_PROTEIN : /* Kimura Proteins*/
							if (!blnNucleo)
								p_KimuraProtein(pmatSites, pmatDistances, pintNbreEspeces, pintNbreSites);
							else
							{	printf("This method is not applicable with nucleotides");
								p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
							return(-1);	}
							break;
	}
	if(ecrireResultats == 1)
		p_EnregistrerResultats(pmatNomsEspeces,pmatDistances,pmatSites,pintNbreEspeces,pintNbreSites,outFile);
	
	return intRetour;
}
// the end !!!

