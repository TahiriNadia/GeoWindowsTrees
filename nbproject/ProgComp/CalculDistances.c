// CalculDistances.cpp: implementation of the CalculDistances class.
//
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>   
#include <ctype.h>
#include <string.h>
#include "global.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

int FALSE=0,TRUE=1;

/****************************************************************************************
   Cette classe est utilisé pour calculer les différentes distances entre plusieurs
   espèces en utilisant l'une des méthodes suivantes :
		1- Uncorrected distances			5- Tamura
		2- Jukes-Cantor						6- Jin-Nei Gamma
		3- Tajima-Nei						7- Kimura Protein
		4- Kimura 2-Parameter				8- LogDet
		9- F84
   
   Il y a des méthodes pour lesquelles le programme demande de l'utilisateur d'entrer 
   des paramètres, et après il lui demande d'entrer un nom de fichier de sortie pour 
   enregistrer les résultats de calcul.

   Remarque intéressante : pour les indices de tableaux et de matrices, on commence par 1 
   au lieu de 0. Donc on ajoute toujours une case supplémentaire lors de l'allocation 
   de mémoire.

   Aussi, on utilise les sites manquants '?' comme des gaps '-'.
*****************************************************************************************/


/****************************************************************************
/*	Nom			:	p_LibererMemoire
/*	Paramètres	:	pmatNomsEspeces (in) : Matrice d'espèces
/*					pmatDistances   (in) : Matrice de distances
/*					pmatSites       (in) : Matrice de sites
/*					pintNbreSites   (in) : Nombre de sites
/*	Retour		:	-
/*	Objectif	:	Libération de la mémoire.
*****************************************************************************/
void p_LibererMemoire(char **pmatNomsEspeces, double **pmatDistances,char**pmatSites, int pintNbreEspeces, int pintNbreSites)
{
	int i;
	if (pintNbreSites != 0)
	{	for (i = 0; i <= pintNbreEspeces; i++)
		{	free(pmatSites[i]); 
			free(pmatDistances[i]); 
		}
		free(pmatSites); free(pmatDistances); free(pmatNomsEspeces);
	}
	else 
	{	for (i = 0; i <= pintNbreEspeces; i++)
			free(pmatDistances[i]); 
		free(pmatDistances); free(pmatNomsEspeces);
	}
}

/***************************************************************************************
/*	Nom			:	f_LireFichierSource
/*	Paramètres	:	pfilPhylip      (in)     : Fichier format Phylip
/*					pmatNomsEspeces (in/out) : Matrice d'espèces
/*					pmatDistances   (in/out) : Matrice de distances
/*					pmatSites       (in/out) : Matrice de sites
/*					pintNbreEspeces (in)     : Nombre d'espèces
/*					pintNbreSites   (in)     : Nombre de sites
/*					pintTailleLigne (in/out) : Taille d'une ligne en terme de caractères
/*					pblnNucleo      (in/out) : Type de données du fichier source (nucléotides/protéines)
/*					pintIdMethod	(in)	 : Identifiant du modèle d'évolution de séquence
/*					pblnWarningDejaAffiche (in) : Indique si on a déjà affiché le warning de Uncorrected distances
/*	Retour		:	0 (OK) ou -1 (KO)
/*	Objectif	:	Lire le fichier source, format Phylip, et remplir la
/*					matrice d'espèces et la matrice de distances et/ou de sites.
*****************************************************************************************/
int f_LireFichierSource(FILE *pfilPhylip, 
                        char **pmatNomsEspeces, 
						double **pmatDistances, 
						char **pmatSites, 
						int pintNbreEspeces, 
						int pintNbreSites, 
						int *pblnNucleo, 
						int pintIdMethod, 
						int pblnWarningDejaAffiche)
{
	
	char chCaratere;		// Caractère lu du fichier source
	int i, j, k;			// Compteurs
	int blnFinLigne;		// Indique le fin d'une ligne
	int blnT = FALSE;		// Indique l'existance du caractère T
	int blnU = FALSE;		// Indique l'existance du caractère U

	char strMsgErreur[200];	// Message d'erruer
	
	int blnPremierProt = FALSE;	// Indiquer le premier caractère de protèine
	char chPremierProt;		// Récupérer le premier caractère de protèine
	int blnWarning = FALSE;	// Indique que le warning lié au modèle uncorrected distances est déjà affiché
	
	if (pintNbreSites == 0) // cas de remplir la matrice des distances
	{	for (i = 1; i <= pintNbreEspeces; i++)
		{	k = 0;
			chCaratere = fgetc(pfilPhylip);
			// Récupérer le nom d'espèce
			while ((chCaratere != ' ') && (chCaratere != '\t'))
			{	if (chCaratere != '\n') 
				{	pmatNomsEspeces[i][k] = chCaratere;
					k++;
				}
				chCaratere = fgetc(pfilPhylip);
			}
			pmatNomsEspeces[i][k] = '\0';

			// Récupérer pintNbreEspeces de distances contenues sur la même ligne du fichier
			for (j = 1; j <= pintNbreEspeces; j++)
				fscanf(pfilPhylip,"%lf", &pmatDistances[i][j]);
		}
	}
	else // cas de remplir la matrice des sites
	{	chCaratere = fgetc(pfilPhylip); // On saute le retour à la ligne
		for (i = 1; i <= pintNbreEspeces; i++)
		{	k = 0;                           
			chCaratere = fgetc(pfilPhylip);
			// Récupérer le nom d'espèce dont la taille ne dépasse pas 20 caractères
			//while (k < 20)
			while ((k < 30) && (chCaratere != '\t') && (chCaratere != ' '))
			{	
				pmatNomsEspeces[i][k] = chCaratere;
				k++;
				chCaratere = fgetc(pfilPhylip);
			}
			pmatNomsEspeces[i][k] = '\0';
		/*	while(k<10){
			//	printf("%d -> %c\n",k,chCaratere);
				if(chCaratere == ' ')
					chCaratere = '_'; 
				pmatNomsEspeces[i][k] = chCaratere;
				k++;
				chCaratere = fgetc(pfilPhylip);
			}
			
			pmatNomsEspeces[i][k] = '\0';
			j=9;
			while (j>0){
				if(pmatNomsEspeces[i][j] == '_')   
					pmatNomsEspeces[i][j] = '\0';
				else
					j=0; 
				j--;
			} */
			
			//printf("\n %s (%d)",pmatNomsEspeces[i],k);
			//exit(8);
			
			if ((k == 30) && (chCaratere != ' ') && (chCaratere != '\t'))
			{	
				printf("Incorrect Phylip format. Object names and sequences\n must be separated by at least one space.");
				return(-1);
			}

			// Dépasser les espaces ainsi que les tabulations
			//chCaratere = fgetc(pfilPhylip);
			while ((chCaratere == ' ') || (chCaratere == '\t'))
				chCaratere = fgetc(pfilPhylip);
			
			// Récupérer les sites (des caractères) jusqu'à la fin de la ligne
			j = 1;
			//while (chCaratere != '\n')
			while ((chCaratere != '\n') && (!feof(pfilPhylip)))
			{	if (chCaratere != ' ')
				{	if (!strchr("-?ARNDCEQGHILKMFPSTWYVUX", toupper(chCaratere))) // x est rajouté le 9 mars 2005 par alix
					{	sprintf(strMsgErreur,"Character %c is not a protein or nucleotide",chCaratere);

						if (pintIdMethod == ID_SM_UNCOR_DIST)
						{	if  ((!blnWarning) && (!pblnWarningDejaAffiche))
							{	printf("%s",strMsgErreur);
								blnWarning = TRUE;
							}
						}
						else
						{	printf("%s",strMsgErreur);
							return(-1);
						}
					}
					else if (strchr("-?AGCTU",toupper(chCaratere)))
					{	if (toupper(chCaratere) == 'T') blnT = TRUE;
						else if (toupper(chCaratere) == 'U') blnU = TRUE;
					}
					else 
					{	*pblnNucleo = FALSE;
						if (!blnPremierProt)
						{	blnPremierProt = TRUE;
							chPremierProt = toupper(chCaratere);
						}
					}
					
					if (blnT && blnU && *pblnNucleo)
					{	if (pintIdMethod == ID_SM_UNCOR_DIST)
						{	if  ((!blnWarning) && (!pblnWarningDejaAffiche))
							{	printf("The characters U and T are both present in this file !");
								blnWarning = TRUE;
							}
						}
						else
						{	//AfxMessageBox("The characters U and T are introduced into the file of nucleotides", MB_ICONERROR);
							printf("Sequence file contains ambiguous characters U and T !\n T-Rex cannot decide whether this file contains DNA or RNA sequences.");
							return(-1);
						}
					}
					//else if (blnU)
					else if (blnU && !*pblnNucleo)
					{	printf("Sequence file contains ambiguous characters U (not a protein) and %c (not a nucleotide).\nT-Rex cannot decide whether this file contains nucleotides or proteins.", chPremierProt);
						return(-1);
					}

					pmatSites[i][j] = chCaratere;
					j++;
				}
				chCaratere = fgetc(pfilPhylip);
			}
		}

		k = j; // Sauvegarder la postion où on est arrivé dans la matrice des sites

		chCaratere = fgetc(pfilPhylip);
		
		while (k <= pintNbreSites)
		{	j = k; // Sauvegarder la postion où on est arrivé dans la matrice des sites
				   // pour les prochaines espèces
			for (i = 1; i <= pintNbreEspeces; i++)
			{	k = j; // Indiquer où on va poursuivre le stockage pour les différents espèces
				// Dépasser les espaces, les tabulations ainsi que les retours à la ligne
				while ((chCaratere == ' ') || (chCaratere == '\t') || (chCaratere == '\n'))
					chCaratere = fgetc(pfilPhylip);
				
				blnFinLigne = FALSE;
				// Remplir la matrice des sites, pour chaque espèce, tant qu'on est pas arrivé
				// à la fin d'une ligne ou à la fin des sites (en terme de caractères)
				while (!blnFinLigne)
				{	if ((chCaratere == '\n') || (k == pintNbreSites + 1)) blnFinLigne = TRUE;
					
					if ((chCaratere != ' ') && (chCaratere != '\t') && (chCaratere != '\n') && (k != pintNbreSites + 1))
					{	if (!strchr("-?ARNDCEQGHILKMFPSTWYVU", toupper(chCaratere)))
						{	sprintf(strMsgErreur,"Character %c is not a protein or nucleotide",chCaratere);
							
							if (pintIdMethod == ID_SM_UNCOR_DIST)
							{	if  ((!blnWarning) && (!pblnWarningDejaAffiche))
								{	printf("%s",strMsgErreur);
									blnWarning = TRUE;
								}
							}
							else
							{	printf("%s",strMsgErreur);
								return(-1);
							}
						}
						else if (strchr("-?AGCTU",toupper(chCaratere)))
						{	if (toupper(chCaratere) == 'T') blnT = TRUE;
							else if (toupper(chCaratere) == 'U') blnU = TRUE;
						}
						else 
						{	*pblnNucleo = FALSE;
							if (!blnPremierProt)
							{	blnPremierProt = TRUE;
								chPremierProt = toupper(chCaratere);
							}
						}
						if (blnT && blnU && *pblnNucleo)
						{	if (pintIdMethod == ID_SM_UNCOR_DIST)
							{	if  ((!blnWarning) && (!pblnWarningDejaAffiche))
								{	printf("The characters U and T are both present in this file !");
									blnWarning = TRUE;
								}
							}
							else
							{
								printf("Sequence file contains ambiguous characters U and T ! T-Rex cannot decide whether this file contains DNA or RNA sequences.");
								return(-1);
							}
						}
						//else if (blnU)
						else if (blnU && !*pblnNucleo)
						{	printf("Sequence file contains ambiguous characters U (not a protein) and %c (not a nucleotide).\nT-Rex cannot decide whether this file contains nucleotides or proteins.", chPremierProt);
							return(-1);
						}

						pmatSites[i][k] = chCaratere;
						k++;
					}
					
					chCaratere = fgetc(pfilPhylip);
				}
			}
		}
	}

	// Mazouzi - début
	/*FILE *dataMazouzi;
	dataMazouzi = fopen("NomsSites.txt","w");
	for (i = 1; i <= pintNbreEspeces; i++)
		fprintf(dataMazouzi, "%s\n", pmatNomsEspeces[i]);
	for (i = 1; i <= pintNbreEspeces; i++)
	{	for (j = 1; j <= pintNbreSites; j++)
			fprintf(dataMazouzi, "%c", pmatSites[i][j]);
		fprintf(dataMazouzi, "\n");
	}
	fclose (dataMazouzi);*/
	// Mazouzi - fin

	return(0);
}

/********************************************************************************
/*	Nom			:	p_UncorrectedDistancesNucleo
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*					pdblPenalite    (in) : Pénalité donnée par l'utilisateur
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Uncorrected distances pour nucléotides, et stoker les distances  
/*					entre elles dans la matrice des distances.
**********************************************************************************/
void p_UncorrectedDistancesNucleo(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
							      int pintNbreSites, double pdblPenalite)
{ 
	int i, j, k;			// Compteurs
	int intCorrespandance;	// Nombre de correspondances entre les sites
	int intPosComparable;	// Nombre de positions comparables, là où il n'y a pas de gaps
	int intGaps;			// Les lacunes contenues dans les séquences

	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	intCorrespandance = 0; 
			intGaps = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] == pmatSites[j][k]) && (pmatSites[i][k] != '-') && (pmatSites[i][k] != '?')) 
					intCorrespandance++;
				else if ((pmatSites[i][k] == '-') || (pmatSites[i][k] == '?') || 
					     (pmatSites[j][k] == '-') || (pmatSites[j][k] == '?'))
					intGaps++;
			}
			intPosComparable = pintNbreSites - intGaps ;
	 
			pmatDistances[i][j] = 1 - intCorrespandance / (intPosComparable + 
														   pdblPenalite * intGaps);
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/********************************************************************************
/*	Nom			:	p_UncorrectedDistancesProt
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*					pdblPenalite    (in) : Pénalité donnée par l'utilisateur
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Uncorrected distances pour protéines, et stoker les distances 
/*					entre elles dans la matrice des distances.
**********************************************************************************/
void p_UncorrectedDistancesProt(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
							int pintNbreSites, double pdblPenalite)
{ 
	int i, j, k;			// Compteurs
	int intCorrespandance;	// Nombre de correspondances entre les sites
	int intPosComparable;	// Nombre de positions comparables, là où il n'y a pas de gaps
	int intGaps;			// Les lacunes contenues dans les séquences

	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	intCorrespandance = 0; 
			intGaps = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] == pmatSites[j][k]) && (pmatSites[i][k] != '-') && (pmatSites[i][k] != '?')) 
					intCorrespandance++;
				else if ((pmatSites[i][k] == '-') || (pmatSites[i][k] == '?') || 
					     (pmatSites[j][k] == '-') || (pmatSites[j][k] == '?'))
					intGaps++;
			}
			intPosComparable = pintNbreSites - intGaps ;
	 
			pmatDistances[i][j] = 1 - intCorrespandance / (intPosComparable + 
														   pdblPenalite * intGaps);
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/*************************************************************************************
/*	Nom			:	p_JukesCantorNucleo
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*					pdblPenalite    (in) : Pénalité donnée par l'utilisateur
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Jukes-Cantor pour nucléotides, et stoker les distances entre elles 
/*					dans la matrice des distances.
**************************************************************************************/
void p_JukesCantorNucleo(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
				         int pintNbreSites, double pdblPenalite)
{ 
	int i, j, k;			// Compteurs
	int intCorrespandance;	// Nombre de correspondances entre les sites
	int intPosComparable;	// Nombre de positions comparables, là où il n'y a pas de gaps
	int intGaps;			// Les lacunes contenues dans les séquences
	double dblD;			// Uncorrected distance
	double dblb;			// Paramètre b

	dblb = 3./4;

	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	intCorrespandance = 0; 
			intGaps = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] == pmatSites[j][k]) && (pmatSites[i][k] != '-') && (pmatSites[i][k] != '?')) 
					intCorrespandance++;
				else if ((pmatSites[i][k] == '-') || (pmatSites[i][k] == '?') ||
						 (pmatSites[j][k] == '-') || (pmatSites[j][k] == '?'))
					intGaps++;
			}
			intPosComparable = pintNbreSites - intGaps ;
			
			dblD = 1 - intCorrespandance / (intPosComparable + pdblPenalite * intGaps);
			dblD = dblD / dblb;
			if (1 - dblD <= 0) pmatDistances[i][j] = '-';
			else pmatDistances[i][j] = - dblb * log(1 - dblD);
			
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/*************************************************************************************
/*	Nom			:	p_JukesCantorProt
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*					pdblPenalite    (in) : Pénalité donnée par l'utilisateur
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Jukes-Cantor pour protéines, et stoker les distances entre elles 
/*					dans la matrice des distances.
**************************************************************************************/
void p_JukesCantorProt(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
				       int pintNbreSites, double pdblPenalite)
{ 
	int i, j, k;			// Compteurs
	int intCorrespandance;	// Nombre de correspondances entre les sites
	int intPosComparable;	// Nombre de positions comparables, là où il n'y a pas de gaps
	int intGaps;			// Les lacunes contenues dans les séquences
	double dblD;			// Uncorrected distance
	double dblb;			// Paramètre b

	dblb = 19./20;

	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	intCorrespandance = 0; 
			intGaps = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] == pmatSites[j][k]) && (pmatSites[i][k] != '-') && (pmatSites[i][k] != '?')) 
					intCorrespandance++;
				else if ((pmatSites[i][k] == '-') || (pmatSites[i][k] == '?') ||
						 (pmatSites[j][k] == '-') || (pmatSites[j][k] == '?'))
					intGaps++;
			}
			intPosComparable = pintNbreSites - intGaps ;
			
			dblD = 1 - intCorrespandance / (intPosComparable + pdblPenalite * intGaps);
			dblD = dblD / dblb;
			if (1 - dblD <= 0) pmatDistances[i][j] = '-';
			else pmatDistances[i][j] = - dblb * log(1 - dblD);
			
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/*************************************************************************************
/*	Nom			:	p_TajimaNei
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Tajima-Nei, et stoker les distances entre elles dans la matrice 
/*					des distances.
**************************************************************************************/
void p_TajimaNei(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
				 int pintNbreSites)
{ 
	int i, j, k, l;				// Compteurs
	double dblCorrespandance;	// Nombre de correspondances entre les sites
	int intPosComparable;		// Nombre de positions comparables, là où il n'y a pas de gaps
	int intGaps;				// Les lacunes contenues dans les séquences
	double dblD;				// Uncorrected distance
	double tabFraction[5];		// Tableau contenant le nombre de substitutions
	double matFrequence[5][5];	// Matrice contenant le nombre de substitutions
	double dblTotalFract;		// Total des fractions
	double dblH;				// Calculé à partir de tabFraction et matFrequence
	double dblB;				// b-parameter


	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	dblCorrespandance = 0; intGaps = 0; dblTotalFract = 0; dblH = 0;
			for (k = 1; k <= 4; k++)
				tabFraction[k] = 0;
			for (k = 1; k <= 3; k++)
				for (l = k+1; l <= 4; l++)
					matFrequence[k][l] = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] == pmatSites[j][k]) && (pmatSites[i][k] != '-') && (pmatSites[i][k] != '?')) 
					dblCorrespandance++;
				else if ((pmatSites[i][k] == '-') || (pmatSites[i][k] == '?') ||
						 (pmatSites[j][k] == '-') || (pmatSites[j][k] == '?'))
					intGaps++;

				// Récupération des fréquences
				if ((pmatSites[i][k] != '-') && (pmatSites[i][k] != '?') &&
					(pmatSites[j][k] != '-') && (pmatSites[j][k] != '?'))
				{	if ((toupper(pmatSites[i][k]) == 'A') || (toupper(pmatSites[j][k]) == 'A'))
					{	tabFraction[1] += 1;
						if (pmatSites[i][k] == pmatSites[j][k]) tabFraction[1] += 1;
					}
					if ((toupper(pmatSites[i][k]) == 'T') || (toupper(pmatSites[j][k]) == 'T'))
					{	tabFraction[2] += 1;
						if (pmatSites[i][k] == pmatSites[j][k]) tabFraction[2] += 1;
					}
					if ((toupper(pmatSites[i][k]) == 'C') || (toupper(pmatSites[j][k]) == 'C'))
					{	tabFraction[3] += 1;
						if (pmatSites[i][k] == pmatSites[j][k]) tabFraction[3] += 1;
					}
					if ((toupper(pmatSites[i][k]) == 'G') || (toupper(pmatSites[j][k]) == 'G'))
					{	tabFraction[4] += 1;
						if (pmatSites[i][k] == pmatSites[j][k]) tabFraction[4] += 1;
					}
				}
				// Récupération des pair
				if (((toupper(pmatSites[i][k]) == 'A') && (toupper(pmatSites[j][k]) == 'T')) ||
					((toupper(pmatSites[i][k]) == 'T') && (toupper(pmatSites[j][k]) == 'A')))
					matFrequence[1][2] += 1;
				else if (((toupper(pmatSites[i][k]) == 'A') && (toupper(pmatSites[j][k]) == 'C')) ||
						 ((toupper(pmatSites[i][k]) == 'C') && (toupper(pmatSites[j][k]) == 'A')))
					matFrequence[1][3] += 1;
				else if (((toupper(pmatSites[i][k]) == 'A') && (toupper(pmatSites[j][k]) == 'G')) ||
						 ((toupper(pmatSites[i][k]) == 'G') && (toupper(pmatSites[j][k]) == 'A')))
					matFrequence[1][4] += 1;
				else if (((toupper(pmatSites[i][k]) == 'T') && (toupper(pmatSites[j][k]) == 'C')) ||
						 ((toupper(pmatSites[i][k]) == 'C') && (toupper(pmatSites[j][k]) == 'T')))
					matFrequence[2][3] += 1;
				else if (((toupper(pmatSites[i][k]) == 'T') && (toupper(pmatSites[j][k]) == 'G')) ||
						 ((toupper(pmatSites[i][k]) == 'G') && (toupper(pmatSites[j][k]) == 'T')))
					matFrequence[2][4] += 1;
				else if (((toupper(pmatSites[i][k]) == 'C') && (toupper(pmatSites[j][k]) == 'G')) ||
						 ((toupper(pmatSites[i][k]) == 'G') && (toupper(pmatSites[j][k]) == 'C')))
					matFrequence[3][4] += 1;
			}
			
			intPosComparable = pintNbreSites - intGaps ;

			tabFraction[1] = tabFraction[1] / (2*intPosComparable);
			tabFraction[2] = tabFraction[2] / (2*intPosComparable);
			tabFraction[3] = tabFraction[3] / (2*intPosComparable);
			tabFraction[4] = tabFraction[4] / (2*intPosComparable);
			
			// Calcul de la moyenne de chaque pair AT/TA AC/CA AG/GA TC/CT TG/GT CG/GC
			for (k = 1; k <= 3; k++)
				for (l = k+1; l <= 4; l++)
					matFrequence[k][l] = matFrequence[k][l] / intPosComparable;
			
			for (k = 1; k <= 3; k++)
				for (l = k+1; l <= 4; l++)
					dblH += (matFrequence[k][l] * matFrequence[k][l]) /
									(tabFraction[k] * tabFraction[l]) ;
			dblH = dblH / 2 ;
			
			for (k = 1; k <= 4; k++)
				dblTotalFract += tabFraction[k] * tabFraction[k];

			dblD = 1 - dblCorrespandance / intPosComparable; // Gaps ignorés

			dblB = (1 - dblTotalFract + (dblD*dblD)/dblH) / 2;
			if (1 - dblD/dblB <= 0) pmatDistances[i][j] = '-';
			else pmatDistances[i][j] = - dblB * log(1 - dblD/dblB);
			
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/*************************************************************************************
/*	Nom			:	p_Kimura2Parameter
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Kimura 2-Parameter, et stoker les distances entre elles dans la  
/*					matrice des distances.
**************************************************************************************/
void p_Kimura2Parameter(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
				        int pintNbreSites)
{ 
	int i, j, k;			// Compteurs
	double dblPosComparable;// Nombre de positions comparables, là où il n'y a pas de gaps
	int intTransition;		// Une transition : A/G, G/A, C/T, T/C, C/U, U/C
	int intTransversion;	// Une transversion : A/C, C/A, A/T, T/A, A/U, U/A
							//					  G/C, C/G, G/T, T/G, G/U, U/G
	double dblP;			// Transitions/positions comparables
	double dblQ;			// Transversions/positions comparables
	double dblD;			// Distance

	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	dblPosComparable = 0; 
			intTransition = 0;
			intTransversion = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] != '-') && (pmatSites[i][k] != '?') &&
					(pmatSites[j][k] != '-') && (pmatSites[j][k] != '?'))
				{	dblPosComparable++;
					if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'G')) || 
						((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'A')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'T')) || 
						((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'C')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'U')) ||  
						((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'C')))
						intTransition++;
					else if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'G')))
						intTransversion++;
				}
			}

			dblP = intTransition / dblPosComparable;
			dblQ = intTransversion / dblPosComparable;
			dblD = (1 - 2*dblP - dblQ) * sqrt(1 - 2*dblQ);
			if ((1 - 2*dblQ < 0) || (dblD <= 0))
				pmatDistances[i][j] = '-';
			else pmatDistances[i][j] = - 0.5 * log(dblD);
			
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/*************************************************************************************
/*	Nom			:	p_Tamura
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Tamura, et stoker les distances entre elles dans la matrice 
/*					des distances.
**************************************************************************************/
void p_Tamura(char **pmatSites, double **pmatDistances, int pintNbreEspeces, int pintNbreSites)
{ 
	int i, j, k;			// Compteurs
	double dblPosComparable;// Nombre de positions comparables, là où il n'y a pas de gaps
	int intTransition;		// Une transition : A/G, G/A, C/T, T/C, C/U, U/C
	int intTransversion;	// Une transversion : A/C, C/A, A/T, T/A, A/U, U/A
							//					  G/C, C/G, G/T, T/G, G/U, U/G
	double dblP;			// Transitions/positions comparables
	double dblQ;			// Transversions/positions comparables
	double dblGC1;			// Nombre de G et C dans la première séquence
	double dblGC2;			// Nombre de G et C dans la deuxième séquence
	double dblC;			// Calculé à partir de intGC1 et intGC2
	
	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	dblPosComparable = 0; intTransition = 0; intTransversion = 0;
			dblGC1 = 0; dblGC2 = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] != '-') && (pmatSites[i][k] != '?') &&
					(pmatSites[j][k] != '-') && (pmatSites[j][k] != '?'))
				{	dblPosComparable++;
					
					if ((pmatSites[i][k] == 'G')  || (pmatSites[i][k] == 'C'))
						dblGC1++;
					if ((pmatSites[j][k] == 'G')  || (pmatSites[j][k] == 'C'))
						dblGC2++;

					if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'G')) || 
						((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'A')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'T')) || 
						((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'C')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'U')) ||  
						((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'C')))
						intTransition++;
					else if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'G')))
						intTransversion++;
				}
			}

			dblP = intTransition / dblPosComparable;
			dblQ = intTransversion / dblPosComparable;

			dblC = (dblGC1 + dblGC2) / dblPosComparable - 
				   2 * (dblGC1 * dblGC2) / (dblPosComparable * dblPosComparable);

			if ((1 - dblP/dblC - dblQ < 0) || (1 - 2*dblQ < 0))
				pmatDistances[i][j] = '-';
			else 
				pmatDistances[i][j] = - dblC * log(1 - dblP/dblC - dblQ) 
				                      - 0.5 * (1 - dblC) * log(1 - 2*dblQ);
										
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/*************************************************************************************
/*	Nom			:	f_JinNeiGamma
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*					pdbla			(in) : Paramètre 'a'
/*	Retour		:	0 (OK) ou -1 (KO)
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Jin-Nei Gamma, et stoker les distances entre elles dans la  
/*					matrice des distances.
**************************************************************************************/
int f_JinNeiGamma(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
				   int pintNbreSites, double pdbla)
{ 
	int i, j, k, l = 1;		// Compteurs
	double dblPosComparable;// Nombre de positions comparables, là où il n'y a pas de gaps
	int intTransition;		// Une transition : A/G, G/A, C/T, T/C, C/U, U/C
	int intTransversion;	// Une transversion : A/C, C/A, A/T, T/A, A/U, U/A
							//					  G/C, C/G, G/T, T/G, G/U, U/G
	int intTransTransv;		// Soit transition, transversion ou aucune des deux
	double dblP;			// Transitions/positions comparables
	double dblQ;			// Transversions/positions comparables
	double **tabMoyenneSub;	// Tableau contenant les moyennes de substitutions
	double dblVariance = 0;	// Variance de la moyenne des substitutions
	double dblA;			// Calculé à partir de tabMoyenneSub et dblVariance
	double dblMoyenneL = 0;	// Moyenne de L
	double dblDistance;		// Distance

	tabMoyenneSub = (double**)malloc((pintNbreEspeces + 1)*sizeof(double*));
	for (i = 0; i <= pintNbreEspeces; i++)
	{	tabMoyenneSub[i] = (double*)malloc((pintNbreEspeces + 1)*sizeof(double));
		/*if (tabMoyenneSub[i] == NULL)
		{	printf("Insufficient memory capacity.");
			return(-1); }*/
	}
	
	if (pdbla == -1)  // On calcule le paramètre 'a'
	{ for (i = 1; i <= pintNbreEspeces - 1; i++)
	  {	for (j = i+1; j <= pintNbreEspeces; j++)
		{	dblPosComparable = 0; intTransition = 0; intTransversion = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] != '-') && (pmatSites[i][k] != '?') &&
					(pmatSites[j][k] != '-') && (pmatSites[j][k] != '?'))
				{	dblPosComparable++;
					if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'G')) || 
						((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'A')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'T')) || 
						((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'C')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'U')) ||  
						((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'C')))
						intTransition++;
					else if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'G')))
						intTransversion++;
				}
			}
			tabMoyenneSub[i][j] = (intTransition + 2*intTransversion) / dblPosComparable;
		}
	  }
	}
	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	dblPosComparable = 0; intTransition = 0; intTransversion = 0;
			dblVariance = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] != '-') && (pmatSites[i][k] != '?') &&
					(pmatSites[j][k] != '-') && (pmatSites[j][k] != '?'))
				{	dblPosComparable++;
					intTransTransv = 0;
					if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'G')) || 
						((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'A')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'T')) || 
						((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'C')) || 
						((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'U')) ||  
						((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'C')))
						{	intTransition++;
							intTransTransv = 1;
						}
					else if (((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'A')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'A')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'C')) || 
						     ((toupper(pmatSites[i][k]) == 'C')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'T')) ||  
							 ((toupper(pmatSites[i][k]) == 'T')  && (toupper(pmatSites[j][k]) == 'G')) || 
							 ((toupper(pmatSites[i][k]) == 'G')  && (toupper(pmatSites[j][k]) == 'U')) || 
						     ((toupper(pmatSites[i][k]) == 'U')  && (toupper(pmatSites[j][k]) == 'G')))
						{	intTransversion++;
							intTransTransv = 2;
						}
					
					if (pdbla == -1)
					  dblVariance += (tabMoyenneSub[i][j] - intTransTransv) * 
					  			     (tabMoyenneSub[i][j] - intTransTransv);
				}
			}

			if (pdbla == -1)
			{	dblVariance = dblVariance / dblPosComparable;
				dblA = (tabMoyenneSub[i][j] * tabMoyenneSub[i][j]) / dblVariance;
			}
			else
				dblA = pdbla;

			dblP = intTransition / dblPosComparable;
			dblQ = intTransversion / dblPosComparable;

			dblDistance = 0.5 * dblA * (pow(1. - (2*dblP) - dblQ,-1./dblA) +
						               (0.5 * pow(1. - (2*dblQ),-1./dblA)) - 1.5);

			if (dblDistance == HUGE_VAL)
				pmatDistances[i][j] = '-';
			else pmatDistances[i][j] = dblDistance;

			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;

	free(tabMoyenneSub);

	return(0);
}

/*************************************************************************************
/*	Nom			:	p_KimuraProtein
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					Kimura Protein, et stoker les distances entre elles dans la matrice 
/*					des distances.
**************************************************************************************/
void p_KimuraProtein(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
				     int pintNbreSites)
{ 
	int i, j, k;				// Compteurs
	double dblCorrespandance;	// Nombre de correspondances entre les sites
	int intPosComparable;		// Nombre de positions comparables, là où il n'y a pas de gaps
	int intGaps;				// Les lacunes contenues dans les séquences
	double dblD;				// Uncorrected distance

	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	dblCorrespandance = 0; 
			intPosComparable = 0; 
			intGaps = 0;

			for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] == pmatSites[j][k]) && (pmatSites[i][k] != '-') && (pmatSites[i][k] != '?')) 
					dblCorrespandance++;
				else if ((pmatSites[i][k] == '-') || (pmatSites[i][k] == '?') ||
						 (pmatSites[j][k] == '-') || (pmatSites[j][k] == '?'))
					intGaps++;
			}
			intPosComparable = pintNbreSites - intGaps ;

			dblD = 1 - dblCorrespandance / intPosComparable;
			dblD = 1 - dblD - 0.2*dblD * dblD;
			if (dblD <= 0) pmatDistances[i][j] = '-';
			else pmatDistances[i][j] = - log(dblD);
			
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}

/************************************************************************
/*	Nom			:	f_CalculerDeterminant
/*	Paramètres	:	pmatSubst (in) : Matrice des substitutions
/*	Retour		:	Déterminant de la matrice pmatSubst ou -1 (KO)
/*	Objectif	:	Calculer le déterminant de la matrice pmatSubst
*************************************************************************/
double f_CalculerDeterminant(double pmatSubst[5][5])
{
	long i, j, k;			// Compteurs
	double dblTemp;			// Temporaire
	double dblDeterminant;	// Déterminant

	dblDeterminant = 1.0;
	for (i = 1; i <= 4; i++) 
	{	dblDeterminant *= pmatSubst[i][i];
		dblTemp = 1.0 / pmatSubst[i][i];
		pmatSubst[i][i] = 1.0;
		for (j = 1; j <= 4; j++)
			pmatSubst[i][j] *= dblTemp;
		for (j = 1; j <= 4; j++) 
		{	if (j != i) 
			{	dblTemp = pmatSubst[j][i];
				pmatSubst[j][i] = 0.0;
				for (k = 1; k <= 4; k++)
				pmatSubst[j][k] -= dblTemp * pmatSubst[i][k];
			}
		}
	}
	
	if (dblDeterminant <= 0.0)
		return(-1);
	else
		return(log(dblDeterminant));
}

/*************************************************************************************
/*	Nom			:	p_LogDet
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					LogDet, et stoker les distances entre elles dans la matrice 
/*					des distances.
**************************************************************************************/
void p_LogDet(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
			  int pintNbreSites)
{ 
	int i, j, k, l;				// Compteurs
	double matSubst[5][5];		// Matrice contenant le nombre de substitutions
	double tabSites1[5];		// Tableau contenant le nombre de sites de la séquence 1
	double tabSites2[5];		// Tableau contenant le nombre de sites de la séquence 2
	double dblDeterminant;		// Déterminant de la matrice des substitutions matSubst

	// Dans cette procédure, les gaps "-" sont pris, dans le progaramme de Phylip, comme 
	// le premier parmi les sites A, C, G et T/U, est qui est donc A.
	
	// On fait ça au lieu de faire plusieurs tests sur les '?', surtout dans cette procédure
	for (i = 1; i <= pintNbreEspeces; i++)
		for (j = 1; j <= pintNbreSites; j++)
			if (pmatSites[i][j] == '?') pmatSites[i][j] = '-'; // '?' est comme un gap '-'
	
	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	for (k = 1; k <= 4; k++)
			{	tabSites1[k] = 0;
				tabSites2[k] = 0;
				for (l = 1; l <= 4; l++)
					matSubst[k][l] = 0;
			}
			
			for (k = 1; k <= pintNbreSites; k++)
			{	
				if ((toupper(pmatSites[i][k]) == 'A') || (pmatSites[i][k] == '-'))
					tabSites1[1] += 1;
				if (toupper(pmatSites[i][k]) == 'C')
					tabSites1[2] += 1;
				if (toupper(pmatSites[i][k]) == 'G')
					tabSites1[3] += 1;
				if (strchr("TU", toupper(pmatSites[i][k])))
					tabSites1[4] += 1;
				if ((toupper(pmatSites[j][k]) == 'A') || (pmatSites[j][k] == '-'))
					tabSites2[1] += 1;
				if (toupper(pmatSites[j][k]) == 'C')
					tabSites2[2] += 1;
				if (toupper(pmatSites[j][k]) == 'G')
					tabSites2[3] += 1;
				if (strchr("TU", toupper(pmatSites[j][k])))
					tabSites2[4] += 1;
				
				// Récupération des substitutions
				if (((toupper(pmatSites[i][k]) == 'A') && (toupper(pmatSites[j][k]) == 'A')) ||
					((toupper(pmatSites[i][k]) == 'A') && (toupper(pmatSites[j][k]) == '-')) ||
					((toupper(pmatSites[i][k]) == '-') && (toupper(pmatSites[j][k]) == 'A')) ||
					((toupper(pmatSites[i][k]) == '-') && (toupper(pmatSites[j][k]) == '-')))
					matSubst[1][1] += 1;
				else if (((toupper(pmatSites[i][k]) == 'A') && (toupper(pmatSites[j][k]) == 'C')) ||
						 ((toupper(pmatSites[i][k]) == '-') && (toupper(pmatSites[j][k]) == 'C')))
					matSubst[1][2] += 1;
				else if (((toupper(pmatSites[i][k]) == 'A') && (toupper(pmatSites[j][k]) == 'G')) ||
						 ((toupper(pmatSites[i][k]) == '-') && (toupper(pmatSites[j][k]) == 'G')))
					matSubst[1][3] += 1;
				else if (((toupper(pmatSites[i][k]) == 'A') && (strchr("TU", toupper(pmatSites[j][k])))) ||
						 ((toupper(pmatSites[i][k]) == '-') && (strchr("TU", toupper(pmatSites[j][k])))))
					matSubst[1][4] += 1;
				else if (((toupper(pmatSites[i][k]) == 'C') && (toupper(pmatSites[j][k]) == 'A')) ||
						 ((toupper(pmatSites[i][k]) == 'C') && (toupper(pmatSites[j][k]) == '-')))
					matSubst[2][1] += 1;
				else if ((toupper(pmatSites[i][k]) == 'C') && (toupper(pmatSites[j][k]) == 'C'))
					matSubst[2][2] += 1;
				else if ((toupper(pmatSites[i][k]) == 'C') && (toupper(pmatSites[j][k]) == 'G'))
					matSubst[2][3] += 1;
				else if ((toupper(pmatSites[i][k]) == 'C') && (strchr("TU", toupper(pmatSites[j][k]))))
					matSubst[2][4] += 1;
				else if (((toupper(pmatSites[i][k]) == 'G') && (toupper(pmatSites[j][k]) == 'A')) ||
						 ((toupper(pmatSites[i][k]) == 'G') && (toupper(pmatSites[j][k]) == '-')))
					matSubst[3][1] += 1;
				else if ((toupper(pmatSites[i][k]) == 'G') && (toupper(pmatSites[j][k]) == 'C'))
					matSubst[3][2] += 1;
				else if ((toupper(pmatSites[i][k]) == 'G') && (toupper(pmatSites[j][k]) == 'G'))
					matSubst[3][3] += 1;
				else if ((toupper(pmatSites[i][k]) == 'G') && (strchr("TU", toupper(pmatSites[j][k]))))
					matSubst[3][4] += 1;
				else if (((strchr("TU", toupper(pmatSites[i][k]))) && (toupper(pmatSites[j][k]) == 'A')) ||
						 ((strchr("TU", toupper(pmatSites[i][k]))) && (toupper(pmatSites[j][k]) == '-')))
					matSubst[4][1] += 1;
				else if ((strchr("TU", toupper(pmatSites[i][k]))) && (toupper(pmatSites[j][k]) == 'C'))
					matSubst[4][2] += 1;
				else if ((strchr("TU", toupper(pmatSites[i][k]))) && (toupper(pmatSites[j][k]) == 'G'))
					matSubst[4][3] += 1;
				else if ((strchr("TU", toupper(pmatSites[i][k]))) && (toupper(pmatSites[j][k]) == 'T'))
					matSubst[4][4] += 1;
			}

			dblDeterminant = f_CalculerDeterminant(matSubst);

			if (dblDeterminant == -1)
				pmatDistances[i][j] = '-';
			else 
				pmatDistances[i][j] = -0.25 * (dblDeterminant - 0.5 * (log(tabSites1[1]) + log(tabSites1[2])
																     + log(tabSites1[3]) + log(tabSites1[4])
																     + log(tabSites2[1]) + log(tabSites2[2])
																     + log(tabSites2[3]) + log(tabSites2[4])));
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
}



/*************************************************************************************
/*	Nom			:	p_RecupererFrequences
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*					pdblFreqA, pdblFreqC, pdblFreqG, pdblFreqT   (in) : Fréquences
/*					pdblFreqR, pdblFreqY, pdblFreqAR    (in) : Fréquences
/*					pdblFreqCY, pdblFreqGR, pdblFreqTY  (in) : Fréquences
/*					pdblRatioTT, pdblXv, pdblFracChange (in) : Intermiediaires
/*	Retour		:	-
/*	Objectif	:	Récupérer des fréquences et certain paramètres pour le modèle F84.
**************************************************************************************/
void p_RecupererFrequences(char **pmatSites, int pintNbreEspeces, int pintNbreSites, double *pdblFreqA, 
						   double *pdblFreqC, double *pdblFreqG, double *pdblFreqT, double *pdblFreqR, 
						   double *pdblFreqY, double *pdblFreqAR, double *pdblFreqCY, double *pdblFreqGR, 
						   double *pdblFreqTY, double *pdblRatioTT, double *pdblXv, double *pdblFracChange)
{
  	int i, j;					// Compteurs
	double dblAA, dblBB, dblXi; // Temporaires
	int intNbreSitesTotal, intGaps = 0;
	char strMessage[100];
	
	*pdblFreqA = 0;
	*pdblFreqC = 0;
	*pdblFreqG = 0;
	*pdblFreqT = 0;

	for (i = 1; i <= pintNbreEspeces; i++)
	{	for (j = 1; j <= pintNbreSites; j++)
		{	if ((pmatSites[i][j] == '-') || (pmatSites[i][j] == '?'))
				intGaps++;
			else if (toupper(pmatSites[i][j]) == 'A') *pdblFreqA += 1;
			else if (toupper(pmatSites[i][j]) == 'C') *pdblFreqC += 1;
			else if (toupper(pmatSites[i][j]) == 'G') *pdblFreqG += 1;
			else if (toupper(pmatSites[i][j]) == 'T') *pdblFreqT += 1;
		}
	}
	intNbreSitesTotal = pintNbreEspeces*pintNbreSites - intGaps;
	*pdblFreqA = *pdblFreqA / intNbreSitesTotal;
	*pdblFreqC = *pdblFreqC / intNbreSitesTotal;
	*pdblFreqG = *pdblFreqG / intNbreSitesTotal;
	*pdblFreqT = *pdblFreqT / intNbreSitesTotal;
	*pdblFreqR = *pdblFreqA + *pdblFreqG;
	*pdblFreqY = *pdblFreqC + *pdblFreqT;
	*pdblFreqAR = *pdblFreqA / *pdblFreqR;
	*pdblFreqGR = *pdblFreqG / *pdblFreqR;
	*pdblFreqCY = *pdblFreqC / *pdblFreqY;
	*pdblFreqTY = *pdblFreqT / *pdblFreqY;

	dblAA = *pdblRatioTT * (*pdblFreqR) * (*pdblFreqY) - *pdblFreqA * (*pdblFreqG) - *pdblFreqC * (*pdblFreqT);
	dblBB = *pdblFreqA * (*pdblFreqGR) + *pdblFreqC * (*pdblFreqTY);

	dblXi = dblAA / (dblAA + dblBB);
	*pdblXv = 1.0 - dblXi;

	if (dblXi < 0.0) 
	{	printf("This transition/transversion ratio is impossible with these base frequencies");
		dblXi = 0.0;
		*pdblXv = 1.0;
		(*pdblRatioTT) = (*pdblFreqA * (*pdblFreqG) + *pdblFreqC * (*pdblFreqT)) / ((*pdblFreqR)*(*pdblFreqY));
		sprintf(strMessage,"Transition/transversion parameter reset\nso transition/transversion ratio is %lf",*pdblRatioTT);
		printf("%s",strMessage);
	}
	if (*pdblFreqA <= 0.0)
		*pdblFreqA = 0.000001;
	if (*pdblFreqC <= 0.0)
		*pdblFreqC = 0.000001;
	if (*pdblFreqG <= 0.0)
		*pdblFreqG = 0.000001;
	if (*pdblFreqT <= 0.0)
		*pdblFreqT = 0.000001;
	
	*pdblFracChange = (dblXi) * (2 * (*pdblFreqA) * (*pdblFreqGR) + 2 * (*pdblFreqC) * (*pdblFreqTY)) +
					  (*pdblXv) * (1.0 - (*pdblFreqA) * (*pdblFreqA) - (*pdblFreqC) * (*pdblFreqC) - (*pdblFreqG) * 
					  (*pdblFreqG) - (*pdblFreqT) * (*pdblFreqT));
}  

/*************************************************************************************
/*	Nom			:	p_F84
/*	Paramètres	:	pmatSites       (in) : Matrice de sites
/*					pmatDistances   (in) : Matrice de distances
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*	Retour		:	-
/*	Objectif	:	Comparer les séquences deux par deux, en utilisant la méthode
/*					F84, et stoker les distances entre elles dans la  matrice des 
/*					distances.
**************************************************************************************/
void p_F84(char **pmatSites, double **pmatDistances, int pintNbreEspeces, 
		   int pintNbreSites)
{ 
	int i, j, k; // Compteurs
	double dblFreqA, dblFreqC, dblFreqG, dblFreqT;	// Fréquences de A, C, G et T
	double dblFreqTempA, dblFreqTempC, dblFreqTempG, dblFreqTempT; // Fréquences temporaires de A, C, G et T
	double dblFreqR, dblFreqY, dblFreqAR, dblFreqCY, dblFreqGR, dblFreqTY; // Fréquences R, Y, AR, CY, GR et TY
	double dblRatioTT = 2.0; // Ratio transition/transvertion
	double dblFracChange; // Fraction de change
	double dblXv, tt, delta, slope, lz, y1, z1yy, z1xv; // Variables temporaires
	int intA, intC, intG, intT; // Variables temporaires
	long it; // Variable temporaire
	
	double dblSomFreq1; // Somme des fréquences de la séquence 1
	double dblSomFreq2; // Somme des fréquences de la séquence 2

	double *tabProd1, *tabProd2, *tabProd3; // Tableaux contenant des produits

	tabProd1 = (double*)malloc((pintNbreSites + 1)*sizeof(double));
	tabProd2 = (double*)malloc((pintNbreSites + 1)*sizeof(double));
	tabProd3 = (double*)malloc((pintNbreSites + 1)*sizeof(double));
	
	p_RecupererFrequences(pmatSites, pintNbreEspeces, pintNbreSites, &dblFreqA, &dblFreqC, &dblFreqG, 
						  &dblFreqT, &dblFreqR, &dblFreqY, &dblFreqAR, &dblFreqCY, &dblFreqGR, 
						  &dblFreqTY, &dblRatioTT, &dblXv, &dblFracChange);
	
	for (i = 1; i <= pintNbreEspeces - 1; i++)
	{	for (j = i+1; j <= pintNbreEspeces; j++)
		{	for (k = 1; k <= pintNbreSites; k++)
			{	if ((pmatSites[i][k] == '-') || (pmatSites[i][k] == '?'))
				{	dblSomFreq1 = 1;
					dblFreqTempA = dblFreqA; dblFreqTempC = dblFreqC; dblFreqTempG = dblFreqG; dblFreqTempT = dblFreqT;
				}
				else if (toupper(pmatSites[i][k]) == 'A')
				{	dblSomFreq1 = dblFreqA;
					dblFreqTempA = dblFreqA; dblFreqTempC = 0; dblFreqTempG = 0; dblFreqTempT = 0;
				}
				else if (toupper(pmatSites[i][k]) == 'C')
				{	dblSomFreq1 = dblFreqC;
					dblFreqTempA = 0; dblFreqTempC = dblFreqC; dblFreqTempG = 0; dblFreqTempT = 0;
				}
				else if (toupper(pmatSites[i][k]) == 'G')
				{	dblSomFreq1 = dblFreqG;
					dblFreqTempA = 0; dblFreqTempC = 0; dblFreqTempG = dblFreqG; dblFreqTempT = 0;
				}
				else if (toupper(pmatSites[i][k]) == 'T')
				{	dblSomFreq1 = dblFreqT;
					dblFreqTempA = 0; dblFreqTempC = 0; dblFreqTempG = 0; dblFreqTempT = dblFreqT;
				}
				
				if ((pmatSites[j][k] == '-') || (pmatSites[j][k] == '?'))
				{	dblSomFreq2 = 1;
					intA = 1; intC = 1; intG = 1; intT = 1;
				}
				else if (toupper(pmatSites[j][k]) == 'A')
				{	dblSomFreq2 = dblFreqA;
					intA = 1; intC = 0; intG = 0; intT = 0;
				}
				else if (toupper(pmatSites[j][k]) == 'C') 
				{	dblSomFreq2 = dblFreqC;
					intA = 0; intC = 1; intG = 0; intT = 0;
				}
				else if (toupper(pmatSites[j][k]) == 'G')
				{	dblSomFreq2 = dblFreqG;
					intA = 0; intC = 0; intG = 1; intT = 0;
				}
				else if (toupper(pmatSites[j][k]) == 'T') 
				{	dblSomFreq2 = dblFreqT;
					intA = 0; intC = 0; intG = 0; intT = 1;
				}
				tabProd1[k] = dblSomFreq1 * dblSomFreq2;
				tabProd2[k] = (dblFreqTempA + dblFreqTempG) * (intA * dblFreqAR + intG * dblFreqGR) +
							  (dblFreqTempC + dblFreqTempT) * (intC * dblFreqCY + intT * dblFreqTY);
				tabProd3[k] = dblFreqTempA * intA + dblFreqTempC * intC +
							  dblFreqTempG * intG + dblFreqTempT * intT;
			}

			tt = 0.1;
			delta = 0.1;
			it = 1;
			while (it < 100 && fabs(delta) > 0.00002) 
			{	slope = 0.0;
				if (tt > 0.0) 
				{	lz = -tt;
					y1 = 1.0 - exp(dblXv * lz);
					z1yy = exp(dblXv * lz) - exp(lz);
					z1xv = exp(dblXv * lz) * dblXv;
					for (k = 1; k <= pintNbreSites; k++)
						slope += (exp(lz) * (tabProd2[k] - tabProd3[k]) + z1xv * (tabProd1[k] - tabProd2[k])) /
								 (tabProd3[k] * exp(lz) + tabProd2[k] * z1yy + tabProd1[k] * y1);
				}
				if (slope < 0.0)
					delta = fabs(delta) / -2.0;
				else
					delta = fabs(delta);
				tt += delta;
				it++;
			}
			if (delta >= 0.1)
				pmatDistances[i][j] = '-';
			else pmatDistances[i][j] = tt * dblFracChange;
			
			pmatDistances[j][i] = pmatDistances[i][j];
		}
		pmatDistances[i][i] = 0;
	}
	pmatDistances[pintNbreEspeces][pintNbreEspeces] = 0;
	
	free(tabProd1);
    free(tabProd2);
    free(tabProd3);
}

/************************************************************************
/*	Nom			:	p_EnregistrerResultats
/*	Paramètres	:	pmatNomsEspeces (in) : Matrice des noms d'espèces
/*					pmatDistances   (in) : Matrice de distances
/*					pmatSites       (in) : Matrice de sites
/*					pintNbreEspeces (in) : Nombre d'espèces
/*					pintNbreSites   (in) : Nombre de sites
/*					pstrMethode     (in) : Nom de la méthode utilisée
/*	Retour		:	-
/*	Objectif	:	Enregistrer dans un fichier de sortie les différentes 
/*					distances entre les séquences comparées.
*************************************************************************/
void p_EnregistrerResultats(char **pmatNomsEspeces, double **pmatDistances, char **pmatSites,
							int pintNbreEspeces, int pintNbreSites, char * strFichierResultat)
{
	FILE *filSauvegarde;			// Fichier de sauvegarde
	char strTabulation[9];			// Ces tabulations servent pour la mise en forme des
									// données écrites dans le fichier de sauvegarde
	int i, j, k;					// Compteurs
	
	if ((filSauvegarde = (FILE*)fopen(strFichierResultat, "w")) == 0)
	{	p_LibererMemoire(pmatNomsEspeces, pmatDistances, pmatSites, pintNbreEspeces, pintNbreSites);
		return;
	}

	fprintf(filSauvegarde, "%d\n", pintNbreEspeces);

	for (i = 1; i <= pintNbreEspeces; i++)
	{	fprintf(filSauvegarde,"%s",pmatNomsEspeces[i]);
		k = strlen(pmatNomsEspeces[i]);
		if (k<=10)
			for(j=k; j<=10; j++)
				fprintf(filSauvegarde," ");
		fprintf(filSauvegarde," ");

		for (j = 1; j <= pintNbreEspeces; j++)
			if (pmatDistances[i][j] == '-')
				fprintf(filSauvegarde,"-99      ");
			else
				//fprintf(filSauvegarde,"%.6f ",pmatDistances[i][j] * 100);
				fprintf(filSauvegarde,"%.6f ",pmatDistances[i][j]);
		fprintf(filSauvegarde,"\n");
	}
	fclose(filSauvegarde);
}
