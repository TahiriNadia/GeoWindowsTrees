//================================================================================================
//=  HGT-DETECTION v3.3.1
//=  Authors : Alix Boc and Vladimir Makarenkov
//=  Date : November 2009
//=
//=  Description : This program detect horizontal gene transfer (HGT). As input it takes 2 
//=  trees: a species tree and a gene tree. the goal is to transform the species tree
//=  into the gene tree following a transfer scenario. There are 3 criteria : the robinson and
//=  Foulds distance, the least-square criterion and the bipartition distance. We also use the
//=  subtree constraint. With this version we can now perform simulation.
//=
//=	 input   : file with species tree and gene tree in the newick format.
//=			   In case of simulation, the species tree and all the gene trees in the same file in
//=            the phylip format or newick string
//=  output  : a list of HGT and the criteria values for each one.
//=	 options :
//=
//================================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#pragma warning(disable:4996)

#include "structures.h"
#include "utils_tree.cpp"
#include "fonctions.cpp"

#define binaireSpecies 0 
#define binaireGene    1

void traiterSignal(int sig){
	printf("\nMESSAGE : SEGMENTATION FAULT #%d DETECTED",sig);
	printf("\nUse valgrind or gdb to fix the problem");
	printf("\n");
	exit(-1);
}

//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
int main(int nargc,char **argv){

	struct InputTree SpeciesTree;				    //== initial species tree
	struct InputTree SpeciesTreeCurrent;		//== initial species tree
	struct InputTree FirstTree;
	struct InputTree FirstTree2;
	struct InputTree FirstTreeB;
	struct InputTree GeneTree;					    //== initial gene tree
	struct InputTree SpeciesTreeRed;			  //== reduced species tree
	struct InputTree GeneTreeRed;				    //== reduced gene tree
	struct ReduceTrace aMap;					      //== mapping structure between species tree and recuded species tree
	struct InputTree geneTreeSave;
	struct HGT * bestHGTRed = NULL;				  //== list of HGT for the reduced tree
	struct HGT * bestHGT = NULL;				    //== list of HGT for the normale tree
	struct HGT * outHGT = NULL;
	struct HGT * bestHGTmulticheck = NULL;
	int nbHGT_boot;
	int first = 1,k,l;
	int cpt_hgt,i,j,tmp,nbTree=0;
	int bootstrap = 0;
	int multigene = 0;
	int nbHgtFound = 0;
	struct CRITERIA * multicheckTab=NULL;
	struct CRITERIA aCrit;						       //== struture of all the criteria
	struct DescTree *DTSpecies,					       //== structure of submatrices for the species tree
					*DTGene;					       //== structure of submatrices for the gene tree
	struct Parameters param;
	FILE *in,*out;
	int max_hgt,nbHGT;
	int ktSpecies;
    int trivial = 1;
	
	int *speciesLeaves = NULL;
	int RFref;
	int imc;	
	char *mot = (char*)malloc(100);
	int nomorehgt=0;
	initInputTree(&geneTreeSave);


	//== read parameters
	//printf("\nhgt : reading options");
	if(readParameters(&param,argv,nargc)==-1){
		printf("\nhgt : no options specified, see the README file for more details\n");
		exit(-1);
	}

	rand_bootstrap = param.rand_bootstrap;
	
	signal(SIGSEGV,traiterSignal);
	
	//== open the input file
	if((in=fopen(param.inputfile,"r"))==NULL){
		printf("\nhgt : The file %s does not exist",param.inputfile);
		exit(-1);
	}
	if(strcmp(param.speciesroot,"file") == 0){
		if(!file_exists(param.speciesRootfileLeaves) && !file_exists(param.speciesRootfile)){
			printf("\nhgt : The file %s does not exist",param.speciesRootfileLeaves);
			exit(-1);
		}
	}
	if(strcmp(param.generoot,"file") == 0){
		if(!file_exists(param.geneRootfileLeaves) && !file_exists(param.geneRootfile)){
			printf("\nhgt : The file %s does not exist",param.geneRootfileLeaves);
			exit(-1);
		}
	}
	if((in=fopen(param.inputfile,"r"))==NULL){
		printf("\nhgt : Cannot open input file (%s)",param.inputfile);
		exit(-1);
	}
	
	initInputTree(&FirstTree);
	initInputTree(&SpeciesTreeCurrent);

	FILE * output;
	if((output = fopen(param.outputfile,"w+"))==NULL){
		printf("PROBLEME AVEC RESULTS");
		exit(0);
	}

	
//==============================================================================
//============================= LECTURE DES ARBRES =============================
//==============================================================================
	
	//printf("\nhgt : reading the input file");
	tmp = readInputFile(in, param.input/*,&SpeciesTree,&GeneTree*/,param.errorFile);

	if(tmp==-1) {
		printf("\nCannot read input data !!\n");
		exit(-1);
	}
  
	cpt_hgt = 0;

	initInputTree(&SpeciesTree);
	initInputTree(&GeneTree);
	initInputTree(&SpeciesTreeRed);
	initInputTree(&GeneTreeRed);

	//== lecture des matrices ou chaines newick en entree
	if(readInput(SPECIE,param.input,&SpeciesTree) == -1){ printf("\nError in species tree\n"); exit(-1);}
	if(readInput(GENE,param.input,&GeneTree) == -1){ printf("\nError in gene tree\n"); getchar(); exit(-1);}

	TrierMatrices(GeneTree.Input,GeneTree.SpeciesName,SpeciesTree.SpeciesName,SpeciesTree.size);
	
	NJ(SpeciesTree.Input,SpeciesTree.ADD,SpeciesTree.size);
	NJ(GeneTree.Input,GeneTree.ADD,GeneTree.size);

	//== construction des differentes représentation des arbres (adjacence,aretes,longueur,degre)
	CreateSubStructures(&SpeciesTree,1,binaireSpecies);
	CreateSubStructures(&GeneTree,1,binaireGene);

	InitCriteria(&aCrit,SpeciesTree.size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	fprintf(output,"%d<>%lf<>%lf\n",aCrit.RF,aCrit.LS,aCrit.BD);

	fclose(output);
	
	return 0;
}
