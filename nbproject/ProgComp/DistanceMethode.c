#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

const char *nucleotid="ACGT";

const char MISSING = '?';
const char GAP = '-';
const char AA = 0;
const char CC = 1;
const char GG = 2;
const char TT= 3;
const char MM = 4;    // pour Missing  
const int NBNUCLEOTIDE = 4;  // the number of nucleotides ACGT

// les methodes de consideration de données manquantes dans le calcul de la 
// matrice de distance.
const int UD = 0;   // uncorrected distance
const int PAUP = 1; // Methode implantée dans paup
const int VB = 2;  // Vladimir Baniré  (nouvelle methode)
const int VB2 = 3; // Vladimir Baniré  (nouvelle methode) utilisant une matrice de poids
const int JC = 4;   // Jukes Cantor
const int TN = 5;   // Tajima Nei
const int K2P= 6;   // Kimura 2 paramètres.

const double NUCWEIGTHS = 3.0/4.0;   // constant that used in JC for nucleotides constant

// valeur de retour des fonctions
const int GOOD = 0;
const int FAIL = 1;

// Valeur de tree generation
const int MAXI = 55;

//static bool useGap;
//static double gap_penalty;

int nbGap(char* s1, char * s2, int taille)
{
	int gap = 0,i;
	for(i = 1; i <= taille; i++)
		if(s1[i] == GAP || s2[i] == GAP)
			gap++;

	return gap;
}

/*****************************************************
* compute number of matching characters between 2 sequences
*
******************************************************/
double PairWMatchUD(char*s1, char* s2, int taille)
{
	//this distance calculation excludes the missing data
	int m = 0;  // score of matches
	int nPos =0;  // number of positions included in m
	int i;

	for(i = 1; i<= taille; i++)
	{
		if((s1[i] != MISSING) && (s2[i] != MISSING)   && (s1[i] != GAP) && (s2[i] != GAP))
		{
			if(s1[i] == s2[i]) m++;
			nPos++;
		}
	
	}
	return (double)m;

}


/*****************************************************
* compute amount of mismatching characters between 2 sequences
*
******************************************************/
double PairWMisMatchUD(char*s1, char* s2, int taille)
{
	//this distance calculation excludes the missing data
	int m = 0;  // score of matches
	int nPos =0;  // number of positions included in m
	int i;
	
	for(i = 1; i<= taille; i++)
	{
		if((s1[i] != MISSING) && (s2[i] != MISSING)&&(s1[i] != GAP) && (s2[i] != GAP))
		{
			if(s1[i] == s2[i]) m++;
			nPos++;
		}
	
	}
	return (double)nPos - (double)m ;

}
/******************************************************
* compute Matches characters between sequences 
*
*******************************************************/
void Compute_UDMatch(char ** C,double ** DI, int n, int taille)
{
	int i,j;
	
	for(i= 1; i <= n; i++)
		for(j = i; j <= n; j++)
		{
			if(i == j)
			{
				DI[i][j] = 0;
			
			}
			else
			{
				DI[i][j] = DI[j][i] = PairWMatchUD(C[i],C[j],taille);
			}
		}
}

/******************************************************
* compute MisMatches characters between sequences 
*
*******************************************************/
void Compute_UDMisMatch(char ** C,double ** DI, int n, int taille)
{
	int i,j;
	
	for(i= 1; i <= n; i++)
		for(j = i; j <= n; j++)
		{
			if(i == j)
			{
				DI[i][j] = 0;
			
			}
			else
			{
				DI[i][j] = DI[j][i] = PairWMisMatchUD(C[i],C[j],taille);
			
			}
		
		}

}

/*******************************************************
* Apply Evolution model using Jukes Cantor
*
********************************************************/

int Evol_JC(char ** C,double ** DI, int n, int taille)
{
	int i,j;
	
	for(i= 1; i <= n; i++)
		for(j = i+1; j <= n; j++)
		{
			//this 2 lines were put to conturn 3/4 mismatching error from JC
		
			
			if( DI[i][j] < 0.74)
				DI[i][j] = DI[j][i] = -NUCWEIGTHS * log(1.0 - (DI[i][j]/NUCWEIGTHS));
			else 
					DI[i][j] = DI[j][i] = '-';
				
		}

	return GOOD;

}
/*********************************************************
* compute distance using the PAUP missing data method
*
**********************************************************/

//this return a code for a nucleotide ACGT? character.
int lettre1(char s)
{
	int r; 
	switch(toupper(s))
	{
	case 'A': r = AA;break;
	case 'C': r =CC;break;
	case 'G': r =GG;break;
	case 'T': r =TT;break;
	//case MISSING: ;break;
	default: r =MM;break;
	
	}

	return r;

}

/**********************************************************
* compute distance using the new method of vladimir &Baniré
*
**********************************************************/
int Compute_VB(char ** C,double ** D, int n, int taille, double gap_penalty)
{



	double *** C1;  // matrix holding substitution probabilities
	int * L;   // number of comparable sites on a specific character
	int * *M;   // number of comparable sites  between sequences

	int i, j, k,l, z;
	int a ; // common divisor of an alignment
	int retr;

    
	// allowing space for L and M
	//********** to be completed
	L=(int *) malloc((taille+1)*sizeof(int));

	M=(int **) malloc((n+2)*sizeof(int*));
	for(i= 0; i < n+2; i++)
		M[i]=(int *) malloc((n+2)*sizeof(int));

	

		// allowing space for C1
	//********** to be completed
	C1=(double ***) malloc((n+2)*sizeof(double**));
	for(i=0; i < n+2; i++)
		C1[i]=(double **) malloc((taille+1)*sizeof(double*));

	for (i= 0; i < n+2; i++)
		for(j= 0; j <= taille; j++) 
		C1[i][j]=(double *) malloc( NBNUCLEOTIDE *sizeof(double));

	
		// initialiser matrice
		for(i = 0; i <n+2 ; i++)
			for (j = 0; j <= taille; j++)
				for(k= 0; k <4; k++)
					C1[i][j][k]= 0.0;
		
	// computing L
	for( i=1; i <= taille; i++)
	{
		a = 0;
		for(j = 1; j <= n; j++)
		{
			if (C[j][i] != MISSING && C[j][i] != GAP)
				a++;	
		}
		L[i] = a; 
	}

	//computing M

	for( i=1; i <= n; i++)
	{
		for(j = 1; j <= n ; j++)
		{
			a = 0;
			for(k = 1; k <= taille; k++)
			{
				if ((C[i][k] != MISSING)  && (C[j][k] != MISSING) 
					&& (C[i][k] != GAP)  && (C[j][k] != GAP))
					a++;	
			}
			M[i][j] = a; 
		}
	}

	// Mazouzi
	//InitMatrice(D,n,n,0,0.0);

	Compute_UDMatch(C,D,n,taille);    // but only matching case

		// computing C1
	for( i=1; i <= n; i++)
	{
		for(j = 1; j <= taille; j++)
		{
			if (L[j] != 0)
			{
				if(C[i][j] == MISSING)
				{
					double t =0;
					for(k = 0; k < 4; k++)  //{ACGT}
					{
						C1[i][j][k] = 0;

						for (l = 1; l <= n; l++)
						{
							if(l != i)   // not the same sequence
							{
								// 2 cases
							
								if(M[i][l] != 0)
								{
									if(C[l][j] != MISSING  &&  C[l][j] != GAP)
									{
										if( C[l][j]  == nucleotid[k])
											C1[i][j][k] += D[i][l]/ (L[j] * M[i][l]);    //  1/a * x/n
										else
											C1[i][j][k] +=  (1.0 / L[j]) * (1.0 / 3.0) *(1.0 - D[i][l]/M[i][l]) ;  // 1/a * (1- x/n) 1/3 			
									}
								}
							}
						}
						t += C1[i][j][k];
					}
				}
				
			}
		}
	}
	
	// Mazouzi
	//InitMatrice(D,n,n,0,0.0);
	Compute_UDMisMatch(C,D,n,taille);    // but only the number of mismatching

	// let's correct D in placing the result in DI

	for( i=1; i <= n; i++)
		for(j= i; j <=n; j++)
		{
			if(i == j )
				D[i][j] =0;
			else
			{
				a = 0;
				for(k = 1; k <= taille; k++)
				{
						if( (C[i][k] == MISSING) && ((C[j][k] != MISSING) &&(C[j][k] != GAP)))
						{
												
								D[i][j] += (1.0 - C1[i][k][lettre1(C[j][k])])/* * (1.0-P[i][j][lettre1(C[j][k])])*/ ;  // remind Position is 0,1,2,3
								a++;
						
						}
						else if( ((C[i][k] != MISSING)  &&(C[i][k] != GAP)) && (C[j][k] == MISSING))
						{
								D[i][j] += (1.0 - C1[j][k][lettre1(C[i][k])])/* * (1.0- P[i][j][lettre1(C[i][k])])*/ ;  // remind Position is 0,1,2,3
								a++;
						
						}
						else if((C[i][k] != MISSING) && (C[j][k] != MISSING) &&  (C[i][k] != GAP) && (C[j][k] != GAP))
						{
							a++;
						}
						else if ((C[i][k] == MISSING) && (C[j][k] == MISSING))
						{
							double x = C1[i][k][0] * (1.0- C1[j][k][0])  + 
								C1[i][k][1] * (1.0- C1[j][k][1]) + C1[i][k][2] * (1.0- C1[j][k][2])
								+ C1[i][k][3] * (1.0- C1[j][k][3]); 

								D[i][j] +=   x; 
								a++;

						// if all both two characters are missing, we ignore them
						}
				
				}
				
				
				if(a == 0)
				{
					D[j][i] = D[i][j] = 1.0;
					
				}
				else
				{
										
						D[i][j] /= a + (nbGap(C[i],C[j],taille) * gap_penalty);   // obtaining the proportion
						D[j][i] = D[i][j];

				}
			}
			
		}

		
	
 
	// released memory for C1 and L,M and P
		
	for (i= 0; i < n+2; i++)
		for(j= 0; j <= taille; j++) 
		 free(C1[i][j]);

	for(i=0; i < n+2; i++)
		free(C1[i]);
	for(i=0; i <n+2; i++)
		free(M[i]);
	
	free(C1);

	

	free(L);

	free(M);

	retr = 	Evol_JC(C,D,n,taille);    // fit to Jukes cantor

	//int retr = 0;
	/*for (i=n; i >1; i--)
		for(j =n ; j > 1; j--)
			D[i][j] = D[i-1][j-1];
	*/
	return retr;

}

/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/

/*******************************************************
* compute distance using Kimura 2 parameters
* using Uncorrected distance approch
*
********************************************************/

// indicate if the character is a purine
int purine(char s){ 
	if((toupper(s) == 'A') || (toupper(s)== 'G'))
		return 1;
	else
		return 0;
}


// indicate if the character is a pyrimidine
int pyrimidine(char s){ 
	if ((toupper(s) == 'C') || (toupper(s)== 'T') || (toupper(s) == 'U'))
		return 1;
	else
		return 0;
}

int transition(char s1, char s2)
{
	if ( (purine(s1) && purine(s2)) || (pyrimidine(s1) && pyrimidine(s2)))
		return 1;
	else
		return 0;
}


int transversion(char s1, char s2)
{
	if ( (purine(s1) && pyrimidine(s2)) || (pyrimidine(s1) && purine(s2)))
		return 1;
	else
		return 0;
}


/*******************************************************
* compute distance using Kimura 2 parameters
* using VB missing data approach
*
********************************************************/
double PairWTransition(char*s1, char* s2, int taille)
{
	//this distance calculation excludes the missing data
	int m = 0;  // score of matches
	int nPos =0;  // number of positions included in m
        int i;
	
	for(i = 1; i<= taille; i++)
	{
		if((s1[i] != MISSING) && (s2[i] != MISSING) && (s1[i] != GAP) && (s2[i] != GAP))
		{
			if((s1[i] != s2[i]) && transition(s1[i],s2[i])) m++;
			nPos++;
		}
	
	}        
	return (double)m;

}

double PairWTransversion(char*s1, char* s2, int taille)
{
	//this distance calculation excludes the missing data
	int m = 0;  // score of matches
	int nPos =0;  // number of positions included in m
        int i;
	
	for(i = 1; i<= taille; i++)
	{
		if((s1[i] != MISSING) && (s2[i] != MISSING) && (s1[i] != GAP) && (s2[i] != GAP))
		{
			if(transversion(s1[i],s2[i])) m++;
			nPos++;
		}
	
	}
	return (double)m;

}


/**********************************************************
* compute distance using the new method of vladimir &Baniré
*
**********************************************************/


int Compute_VB_K2P(char ** C,double ** D, int n, int taille)
{

	double *** C1;  // matrix holding substitution probabilities
	int * L;   // number of comparable sites on a specific character
	int * *M;   // number of comparable sites  between sequences
	
	int i, j, k,l, z;
	int a ; // common divisor of an alignment

	int retour = GOOD;   // désigne la valeur de retour pour la fonction
        
	double P, Q; 
	
	// allowing space for L and M
	//********** to be completed
	L=(int *) malloc((taille+1)*sizeof(int));

	M=(int **) malloc((n+2)*sizeof(int*));
	for(i= 0; i < n+2; i++)
		M[i]=(int *) malloc((n+2)*sizeof(int));

	

		// allowing space for C1
	//********** to be completed
	C1=(double ***) malloc((n+2)*sizeof(double**));
	for(i=0; i < n+2; i++)
		C1[i]=(double **) malloc((taille+1)*sizeof(double*));

	for (i= 0; i < n+2; i++)
		for(j= 0; j <= taille; j++) 
		C1[i][j]=(double *) malloc( NBNUCLEOTIDE *sizeof(double));


	
		// initialiser matrice
		for(i = 0; i <=n ; i++)
			for (j = 0; j <= taille; j++)
				for(k= 0; k <4; k++)
					C1[i][j][k]= 0.0;
	// computing L
	for( i=1; i <= taille; i++)
	{
		a = 0;
		for(j = 1; j <= n; j++)
		{
			if (C[j][i] != MISSING && C[j][i] != GAP)
				a++;	
		}
		L[i] = a; 
	}

	//computing M

	for( i=1; i <= n; i++)
	{
		for(j = 1; j <= n ; j++)
		{
			a = 0;                            
			for(k = 1; k <= taille; k++)
			{
				if ((C[i][k] != MISSING)  && (C[j][k] != MISSING)
					&& (C[i][k] != GAP)  && (C[j][k] != GAP))
					a++;	
			}
			M[i][j] = a; 
		}
	}

	// Mazouzi
	//InitMatrice(D,n,n,0,0.0);

	
	
	//Compute_UD_K2PMatch(C,D,n,taille);    // but only matching case
	Compute_UDMatch(C,D,n,taille);
	

	// computing C1
	for( i=1; i <= n; i++)
	{
		for(j = 1; j <= taille; j++)
		{
			if (L[j] != 0)
			{
				if(C[i][j] == MISSING)
				{
					double t =0;
					for(k = 0; k < 4; k++)  //{ACGT}
					{
						C1[i][j][k] = 0;

						for (l = 1; l <= n; l++)
						{
							if(l != i)   // not the same sequence
							{
								// 2 cases
							
								if(M[i][l] != 0)
								{
									if(C[l][j] != MISSING && C[l][j] != GAP)
									{
										if( C[l][j]  == nucleotid[k])
											C1[i][j][k] += D[i][l]/ (L[j] *M[i][l]);    //  1/a * x/n
										else
											C1[i][j][k] +=  (1.0 / L[j]) * (1.0 / 3.0) *(1.0 - D[i][l]/M[i][l]) ;  // 1/a * (1- x/n) 1/3 			
									}
								}
							}
						}
						t += C1[i][j][k];
					}
				}
				
			}
		}
	}
	
	// Mazouzi
	//InitMatrice(D,n,n,0,0.0);
	
	// let's correct D in placing the result in DI

	      
	for( i=1; i <= n; i++)
		for(j= i; j <=n; j++)
		{
			if(i == j )
				D[i][j] =0;
			else
			{
				P = PairWTransition(C[i],C[j],taille);
				Q = PairWTransversion(C[i],C[j],taille);

				a = M[i][j];
				for(k = 1; k <= taille; k++)
				{
					if (L[k] != 0) 
					{
						if( (C[i][k] == MISSING) && (C[j][k] != MISSING) && (C[j][k] != GAP))
						{
								
								P+= C1[i][k][(lettre1(C[j][k]) +2) % 4];
								Q+= C1[i][k][(lettre1(C[j][k]) +1) % 4] 
									+ C1[i][k][(lettre1(C[j][k]) +3) % 4];
								a++;


						}
						else if( (C[i][k] != MISSING) && (C[i][k] != GAP) && (C[j][k] == MISSING))
						{
						    
								P+= C1[j][k][(lettre1(C[i][k]) +2) % 4];
								Q+= C1[j][k][(lettre1(C[i][k]) +1) % 4] 
									+ C1[j][k][(lettre1(C[i][k]) +3) % 4];
								a++;

						}
						else if((C[i][k] != MISSING) && (C[j][k] != MISSING) 
							&& (C[i][k] != GAP) && (C[j][k] != GAP))
						{
							
						}
						else if ((C[i][k] == MISSING) && (C[j][k] == MISSING))
						{
							
						// if all both two characters are missing, we ignore them
							P += C1[i][k][0] * C1[j][k][2] + C1[i][k][1] * C1[j][k][3]
							   + C1[i][k][2] * C1[j][k][0] + C1[i][k][3] * C1[j][k][1];

							Q += C1[i][k][0] * (C1[j][k][1] + C1[j][k][3]) 
							   + C1[i][k][1] * (C1[j][k][0] + C1[j][k][2])
							   + C1[i][k][2] * (C1[j][k][1] + C1[j][k][3])
							   + C1[i][k][3] * (C1[j][k][0] + C1[j][k][2]); 
							a++;
							
						}
					}
				
				}
				
				
				if(a == 0)
				{
					D[i][j] = D[j][i] = '-';
				}
				else
				{
				
						 P /=  (double) a;
						Q /= (double) a;

					D[j][i] = D[i][j] =	 -0.5 * log((1.0-2.0*P-Q)*sqrt(1.0-2.0*Q));
						if(D[i][j] < 0 )
					{
						D[j][i] = D[i][j] = '-';
					}
					
	
				}
				
			}
			
		}



	/*for (i=n; i >1; i--)
		for(j =n ; j > 1; j--)
			D[i][j] = D[i-1][j-1];
*/
	// released memory for C1 and L,M and P
		
	for (i= 0; i < n+2; i++)
		for(j= 0; j <= taille; j++) 
		 free(C1[i][j]);

	for(i=0; i < n+2; i++)
		free(C1[i]);
	
	for(i=0; i <n+2; i++)
		free(M[i]);
	

	free(C1);

	

	free(L);

	free(M);

	
	return retour;

}
