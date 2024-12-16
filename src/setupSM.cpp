/* Copyright (C) 2009-2015 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* STRINGMOL is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


/* Description
 *
 * This file contains helper functions to set up stringmol runs
 *
 * */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mt19937-2.h"
#include "randutil.h"
#include "params.h"
//}

//stringmol
#include "alignment.h"
#include "agent.h"
#include "SMspp.h"

//metabolism
#include "rules.h"
#include "opcodes.h"
#include "agents_base.h"
#include "stringPM.h"


// Writing PNGs
#include "lodepng.h"
#include <iostream>


#include "setupSM.h"

/* Used in comass_ga and comass_ga_boostwinners
 *
 */
void clearfiles(char *argv[]){

	char fn[256];
	FILE *ftmp;

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	ftmp = fopen(fn,"w");
	fclose(ftmp);
	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);
}



/* Used in comass_ga and comass_ga_boostwinners
 *
 */
void setupSMol(struct runparams &RunPar, int argc, char *argv[]){

	FILE *fp;

	unsigned int rerr,rlim=20;
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = ParameterReadUnsignedInt(fp,"NTRIALS",&rlim,1);
		switch(rerr){
		case 2:
			printf("Multiple NTRIALS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NTRIALS to %u\n",rlim);
			break;
		default:
			printf("NTRIALS not specified;\nSetting NTRIALS to %u\n",rlim);
			break;
		}
		fclose(fp);
	}

	//Read nsteps:
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = ParameterReadUnsignedInt(fp,"NSTEPS",&(RunPar.maxnsteps),1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %u\n",RunPar.maxnsteps);
			RunPar.indefinite=0; //TODO: fix the indefinite thing if NSTEPS is not specified..
			break;
		default:
			printf("NSTEPS not specified;\nEach trial will run to extinction\n");
			break;
		}
		fclose(fp);
	}
	fflush(stdout);

	return;

}



/* Used in comass_ga and comass_ga_boostwinners
 *
 */
void record_spp(stringPM *A){
	char fn[256];
	FILE *fp;

	A->SpeciesPrintCount(stdout,0,-1);

	//printf("Printing species list\n");
	sprintf(fn,"splist%03d.dat",A->run_number);
	if((fp = fopen(fn,"w"))!=NULL){
		A->spl->SpeciesListPrint(fp);
		fclose(fp);
	}
	else{
		printf("Unable to write to %s\n",fn);
	}
}

/* used in SmPm_conpop
 *
 */
float ctspp(stringPM *A, const int spp){

	float count = 0;

	s_ag *pa;


	for(pa = A->nowhead; pa!=NULL; pa = pa->next){
		if(pa->spp->spp==spp)
			count++;
	}

	return count;

}





/*******************************************************************************
* @brief set up the popdy output file
*
* @param[in] A the stringPM object (i.e. the "bucket")
*
* @param[in] overwrite whether to overwrite or not... TODO(sjh): investigate!
*
* @return a double between 0 and 1
*******************************************************************************/
void PopdyInitFile(stringPM *A, bool overwrite){
	char pfn[128];
	memset(pfn,0,128*sizeof(char));

	FILE *ftmp;

	sprintf(pfn,"popdy%03d.dat",A->run_number);

	if(!overwrite){
		FilenameGetUnused(pfn);
	}

	//Make sure this file is empty...
	ftmp = fopen(pfn,"w");
	fclose(ftmp);

	//Record the file name so we can append to it later
	strcpy(A->popdyfn,pfn);

}





/*******************************************************************************
* @brief print the count of each species in the present timestep
*
* @details used in comass_AlifeXII, energetic_AlifeXII, origlife,
*          comass_GA, comass_GA_boostwinners, SmPm_AlifeXII,
*          SmPm_conpop and speigmonst
*
* @param[in] A the stringPM object (i.e. the "bucket")
*
* @param[in] seed used to seed the rng. if <0, chosen by the program
*
* @return a double between 0 and 1
*******************************************************************************/
void SpeciesPrintCounts(stringPM *A, int timestep){

	//char fn[128];
	FILE *fp;
	s_ag *pa;
	int spc,count;
	int finished = 0;
	int nag,*done;

	nag = A->AgentsCount(A->nowhead,-1);

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	//TODO: need to check that this isn't going to cause problems with callers other than smspatial...!
	fp = fopen(A->popdyfn,"a");

	do{
		int i = 0;
		int found=0;
		finished = 1;
		for(i=0,pa=A->nowhead;i<nag;i++,pa=pa->next){
			if(!done[i]){
				if(!found){
					done[i]=1;
					count=1;
					finished=0;
					found=1;
					spc = pa->spp->spp;
				}
				else{
					if(pa->spp->spp==spc){
						done[i]=1;
						count++;
					}
				}
			}
		}

		//Write to file
		if(!finished)
			fprintf(fp,"%d,%d,%d\n",timestep,spc,count);

	}while(!finished);

	fflush(fp);
	fclose(fp);
	free(done);
}




int run_one_comass_trial(const int rr, stringPM *A,  int * params, struct runparams *R){


	int *maxcode;

	//todo: what is the relationship between maxcode and params...?
	//if(!rr)
	maxcode = (int *) malloc(A->blosum->N * sizeof(int));
	memset(maxcode,0,A->blosum->N*sizeof(int));


	A->run_number=rr;
	PopdyInitFile(A);

	//todo DELETE if we don't need epochs any more
	//int lastepoch=A->get_ecosystem(),thisepoch,nepochs=1;

	A->domut=0;
	unsigned int nsteps=0;
	int i;
	for(i=0;R->indefinite || nsteps <= R->maxnsteps;i++){

		A->timestep = i;

		A->comass_TimestepIncrement();
		A->UpdateNowNext();


		if(!(i%1000)){
			record_spp(A);

			//todo: put the below in a function
			printf("%03d At  time %d e=%d\n",rr,i,(int)A->energy);
			SpeciesPrintCounts(A,i);

			setmaxcode(A,maxcode);

			/*
			printf("CODE\tMAX\tMASS\tMAX_USED\n");
			for(int k=0;k<A->blosum->N;k++){
				printf("%c:\t%d\t%d\t%d",A->blosum->key[k],params[k],A->mass[k],maxcode[k]);
				if(A->mass[k]<0 || A->mass[k]>params[k])
					printf(" ERROR\n");
				else
					printf("\n");

			}
			printf("\n");
			*/
		}

		if(!A->AgentsCount(A->nowhead,-1)){// || (!rr && nsteps >1500000) || nsteps >15000000){
			printf("DEATH\n");
			printf("At  time %d e=%d\t",i,(int)A->energy);
			A->SpeciesPrintCount(stdout,0,-1);
			//nsteps=i;
			break;
		}

		nsteps++;
		A->energy += 20;
	}

	free(maxcode);
	return i;
}


void setmaxcode(stringPM *A, int *maxcode){

	s_ag *pag;


	for(int i=0;i<A->blosum->N;i++){
		int count=0;
		for(pag=A->nowhead;pag!=NULL;pag=pag->next){
			int C = strlen(pag->S);
			for(int c=0; c<C; c++){
				if(pag->S[c]==A->blosum->key[i]){
					count++;
				}
			}
		}
		maxcode[i]=count<maxcode[i]?maxcode[i]:count;
	}
}

//These are the functions we need to manipulate mutation networks:
//int ** random_mtx(const int N){
//	int ** matrix;
//
//
//}


void setmutnet(const int * mutnet, swt *blosum){

	int i,j;
	for(i=0;i<blosum->N;i++)
		for(j=0;j<blosum->N;j++){
			blosum->adj[i][j] = mutnet[(i*blosum->N)+j];

		}
}





/******************************************************************************
* @brief load config specified by command-line arguments
*
* @param[in] A the stringPM object
*
* @param[in] argc from "main()"
*
* @param[in] argv from "main()"
*
* @param[in] verbose verbose output
*
* @return 0 if error; 1 if success... todo reverse
*****************************************************************************/
int ParametersLoadFromMainArgs(stringPM *A, int argc, char *argv[], int verbose ){

	//TODO: This can be a bit fragile for runs with >2 arguments... be careful!
	switch(argc){
	case 3:
		if(verbose)printf("Traditional config\n");
		A->ParametersLoad(argv[2],0,1);
		A->AgentsLoad(argv[2],NULL,0,0);
		return 1;
	case 4:
		if(verbose)printf("ATTENTION! loading with separate .mtx file!\n");
		A->ParametersLoad(argv[2],0,1);
		A->AgentsLoad(argv[2],argv[3],0,0);
		return 1;
	default:
		if(verbose)printf("Error: wrong number of arguments - try 2 or 3\n");
		return 0;
	}


}





void print_params(stringPM *A, int ntrials, int nsteps){

	printf("Non-stringPM variables:\n");
	if(ntrials<0)
		printf("NTRIALS     not set - the default value would be used if needed\n");
	else
		printf("NTRIALS     %d\n",ntrials);

	if(nsteps<0)
		printf("NSTEPS      not set - the default value would be used if needed\n");
	else
		printf("NSTEPS      %d\n",nsteps);

	//load params:
	printf("CELLRAD     %f (vcellrad = %f)\n",A->cellrad,A->vcellrad);
	printf("AGRAD       %f\n",A->agrad);
	printf("ENERGY      %d\n",(int) A->energy);
	printf("NSTEPS      %f\n",A->nsteps);

	//load agents: load_table(_mtx); load_mut; load_decay
	if(A->blosum == NULL)
		printf("BLOSUM      not set - needs to be loaded explicitly\n");
	else
		printf("BLOSUM      %d size table loaded\n",A->blosum->N);
	printf("MUTATE      indelrate = %f; subrate = %f\n",A->indelrate,A->subrate);
	printf("DECAY       %f\n",A->decayrate);
	printf("MAXLEN      %u, (maxl0 = %u)\n",A->maxl, A->maxl0);
	printf("ESTEP       %u\n",A->estep);


}





/* This is used in the comass GA - see test.cpp for how to set up randseed properly
 *
 */
void init_randseed_config(int argc, char *argv[]){

	/*TODO: set up random seed properly - the best result was *without* a random seed, and has been lost */
	/* Funny business with unsigned longs... */
	//printf("%u\n", -873302838);//-2041524348);
	char test[128];
	memset(test,0,128*sizeof(char));
	sprintf(test,"-873302838");
	long dummy;
	sscanf(test,"%ld",&dummy);
	//printf("Random seed is %d\n",dummy);
	//printf("Random seed is %lu (unsigned)\n",dummy);


	unsigned long seedin =  2846144656u;

	unsigned int qnnscoring = 1;

	FILE *fpr;
	if((fpr=fopen(argv[2],"r"))!=NULL){
		unsigned int stmp;
		int rerr = ParameterReadUnsignedInt(fpr,"RANDSEED",&stmp,1);
		//TODO: load the full RNG state using load_mt (RNGFILE in config)
		if(rerr){
			printf("Error reading randseed\n");
		}


		rerr = ParameterReadUnsignedInt(fpr,"GAQNN",&qnnscoring,1);
		if(rerr)
			qnnscoring=1;
		seedin = stmp;
		fclose(fpr);
	}

	unsigned long rseed = RandomInitLong(&seedin);//437);//-1);//437);
	//unsigned long rseed = longinitmyrand(NULL);//437);//-1);//437);
	FILE *frs;
	if((frs=fopen("randseed.txt","w"))==NULL){
		printf("Coundln't open randseed.txt\n");
		getchar();
	}else{
		fprintf(frs,"(unsigned) random seed is %lu \n",rseed);
		fflush(frs);
		fclose(frs);
	}

}





/*
 * NB: To get identical trials to those run for ALifeXII, do the following:
 * 1: Use the file "replicase.conf" as the input
 * 2: Fix the random number seed to 437
 * 3: #define DO_ANCESTRY to get ancestry files out...
 * 4: Run on a 32-bit linux slackware system, circa 2008 vintage...:)
 */
//todo: SmPm_AlifeXII() should call this
int run_one_AlifeXII_trial(stringPM *A){

	int i;

	A->AgentsPrint(stdout,"NOW",0);
	A->run_number=0;

	int nsteps=0;
	//TODO: Accommodate indefinitre runs, like this:
	//for(i=0;indefinite || nsteps <= maxnsteps;i++){

	for(i=0;nsteps <= A->nsteps;i++){

		//TODO: find out what this does - rename the variable to  make it clear.
		A->timestep = i;

		A->TimestepIncrement();
		A->UpdateNowNext();

		if(!(i%1000)){
			A->SpeciesPrintCount(stdout,0,-1);
		//}
		//
		//if(!(i%1000)){
			printf("At  time %d e=%d, mutrate = %0.9f & %0.9f\n",i,(int)A->energy,A->subrate,A->indelrate);
			SpeciesPrintCounts(A,i);
		}

//TODO: See equivalent line in SmPm_AlifeXII() function for what should be in the following #ifdef...
//#ifdef DO_ANCESTRY
//#endif

		if(!A->AgentsCount(A->nowhead,-1)){
			printf("DEATH\n");
			printf("At  time %d e=%d, mutrate = %0.9f & %0.9f\t",i,(int)A->energy,A->indelrate,A->subrate);
			A->SpeciesPrintCount(stdout,0,-1);
			//nsteps=i;
			break;
		}
		nsteps++;

		A->energy += A->estep;
	}

	printf("Finished - alls well!\nclear out memory now:\n");
	fflush(stdout);

	//TODO: need to do this outside the function!
	//A.clearout();

	return 0;
}





/*******************************************************************************
* @brief select a random cell in the moore neighbourhood
*
* @param[in] X the x position of the cell
*
* @param[in] Y the y position of the cell
*
* @param[in] Xlim the X dimension of the grid
*
* @param[in] Ylim the Y dimension of the grid
*
* @param[in] xout pointer to the x position of the neighbour
*
* @param[in] yout pointer to the y position of the neighbour
*
* @return 0 always
*******************************************************************************/
int GridSelectRandomMooreNeighbour(const int X, const int Y, const int Xlim, const int Ylim, int *xout, int *yout){

	int pos = 8. * RandomBetween0And1();

	/*    012
	 *    3*4
	 *    567
	 */

	int xoff = 0;
	int yoff = 0;

	switch(pos){
	case 0:
		xoff = -1;
		yoff = -1;
		break;
	case 1:
		yoff = -1;
		break;
	case 2:
		xoff = 1;
		yoff = -1;
		break;
	case 3:
		xoff = -1;
		break;
	case 4:
		xoff = 1;
		break;
	case 5:
		xoff = -1;
		yoff = 1;
		break;
	case 6:
		yoff = 1;
		break;
	case 7:
		xoff = 1;
		yoff = 1;
		break;
	}

	*xout = (X+Xlim+xoff)%Xlim;
	*yout = (Y+Ylim+yoff)%Ylim;
	return 0;
}





//todo(sjh): check whether x and y are members of A - delete if so!
/*******************************************************************************
* @brief find a partner to react with (if one exists) on a grid
*
* @param[in] A the current agent
*
* @param[in] run the grid data
*
* @param[in] x the x position of the agent
*
* @param[in] y the y position of the agent
*
* @return mt_error_code
*******************************************************************************/
s_ag * ReactionSeekRandomSpatialPartner(stringPM *A, smsprun *run,int x, int y){

	int i,j,xx,yy;

	//first, let's count the agents
	int count = 0;
	for(i=-1;i<2;i++){
		for(j=-1;j<2;j++){
			if( !(i == 0 && j ==0) ){
				xx = (x + i + run->gridx)%run->gridx;
				yy = (y + j + run->gridy)%run->gridy;
				if(run->grid[xx][yy]!=NULL){
					if(run->grid[xx][yy]->status == B_UNBOUND){
						if(run->status[xx][yy] == G_NOW){
							count++;
						}
					}
				}
			}
		}
	}
	if(!count)
		return NULL;

	int it = count * RandomBetween0And1();


	//now, let's choose the agents
	count = 0;
	for(i=-1;i<2;i++){
		for(j=-1;j<2;j++){
			if( !(i == 0 && j ==0) ){
				xx = (x + i + run->gridx)%run->gridx;
				yy = (y + j + run->gridy)%run->gridy;
				if(run->grid[xx][yy]!=NULL){
					if(run->grid[xx][yy]->status == B_UNBOUND){
						if(run->status[xx][yy] == G_NOW){
							if(count==it)
								run->status[xx][yy] = G_NEXT;
							return run->grid[xx][yy];
							//count++;
						}
					}
				}
			}
		}
	}

	//We should never get to here! todo(sjh): exit with error code
	printf("Something's wrong - neighbour detected in first pass but none selected\n");
	return NULL;
}





/*******************************************************************************
* @brief set the 'next' grid entries to 'now'
*
* @param[in] run the grid info
*******************************************************************************/
void TimestepGridIncrement(smsprun *run){
	for(int i=0;i<run->gridx;i++){
		for(int j = 0;j<run->gridy;j++){
			if(run->grid[i][j]!=NULL){
				if(run->status[i][j]==G_NEXT)
					run->status[i][j] = G_NOW;
			}
			else{
				run->status[i][j] = G_EMPTY;
			}
		}
	}
}





/*******************************************************************************
* @brief position an agent on the grid
*
* @param[in] ag the agent
*
* @param[in] x the x position
*
* @param[in] y the y position
*******************************************************************************/
void AgentPlaceOnGrid(s_ag *ag,smsprun *run,int x, int y){
	//TODO: error checking!
	run->grid[x][y] = ag;
	//Set status to next - placing the molecule has used up the 'now'
	run->status[x][y] = G_NEXT;
	ag->set=true;
	ag->x=x;
	ag->y=y;
}





/*******************************************************************************
* @brief place a new agent in the Moore neighbourhood
*
* @param[in] A the stringmol bucket
*
* @param[in] run the grid info
*
* @param[in] c the agent
*
* @param[in] x the x position
*
* @param[in] y the y position
*
* @return c if placed, else NULL
*******************************************************************************/
s_ag * AgentPlaceInMooreNeighbourhood(stringPM *A, smsprun *run,s_ag *c,int x,int y){


	int xx,yy;
	int nvacant=0;
	for(int i=-1;i<2;i++){
		for(int j=-1;j<2;j++){
			xx = (x+i+run->gridx)%run->gridx;
			yy = (y+j+run->gridy)%run->gridy;
			//No need to ingnore x,y because it is occupied by the parent!
			if(run->grid[xx][yy]==NULL)
				nvacant++;
		}
	}
	if(nvacant){
		//Decide where to put it:
		int pos = nvacant * RandomBetween0And1();
		int here=0;
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				xx = (x+i+run->gridx)%run->gridx;
				yy = (y+j+run->gridy)%run->gridy;
				//No need to ingnore x,y because it is occupied by the parent!
				if(run->grid[xx][yy]==NULL){
					if(here == pos){
						AgentPlaceOnGrid(c,run,xx,yy);
						return c;
					}
					here++;
				}
			}
		}
	}
	else{
		return NULL;
	}
	return NULL;
	/*TODO(sjh): Need to decide whether to replace or not
	else{
		int noccupied = 8-nvacant;
		if(run->grid[x][y]==NULL){
			//This should never happen....
			printf("Alert! empty parent cell!\n");
			noccupied++;
		}

		int pos = nvacant * rand0to1();
		int here=0;
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				xx = (x+i+run->gridx)%run->gridx;
				yy = (y+j+run->gridy)%run->gridy;
				//No need to ingnore x,y because it is occupied by the parent!
				if(i!=0 && j!=0){
					if(run->grid[xx][yy]!=NULL){
						if(here == pos){
							//Remove the incumbent


							run->grid[xx][yy]=c;
							return;
						}
						here++;
					}
				}
			}
		}
	}
	*/
}





/*******************************************************************************
* @brief cleave
*
* @param[in] A the stringmol bucket
*
* @param[in] run the grid info
*
* @param[in] act the agent
*
* @return which agents have been destroyed (if any) 0: none; 1: active only
 *         2: passive only; 3: both
*******************************************************************************/
int OpcodeCleaveSpatial(stringPM *A, smsprun *run, s_ag *act){//, int x, int y){

	int dac = 0,cpy;
	s_ag *c,*pass,*csite;
	c = NULL;
	pass = act->pass;

	//pick the mol containing the cleave site:
	csite = act->ft?act:pass;

	if(act->f[act->ft]-csite->S < csite->len){

		//1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT
		c = A->AgentMake(pass->label);//,1);

		//Copy the cleaved string to the agent
		char *cs;
		c->S =(char *) malloc(A->maxl0*sizeof(char));
		memset(c->S,0,A->maxl0*sizeof(char));
		cs = csite->S;
		cpy = strlen(cs);

		//Check that we aren't creating a zero-length molecule:
		if(!cpy){
			printf("WARNING: Zero length molecule being created!\n");
		}

		cpy -= act->f[act->ft]-cs;

		if(!cpy){
			printf("ERROR: Zero length molecule definitely being created!\nbail..\n");
			A->AgentFree(c);
			c=NULL;
		}
		else{

			strncpy(c->S,act->f[act->ft],cpy);
			c->len = strlen(c->S);

#ifdef VERBOSE
		printf("String %d created:\n%s\n",c->idx,c->S);
#endif

			//Check the lineage
			A->SpeciesListUpdate(c,'C',1,act->spp,pass->spp,act->biomass);
			act->biomass=0; //reset this; we might continue to make stuff!

			//TODO: place the new agent on the grid
			if((AgentPlaceInMooreNeighbourhood(A,run,c,act->x,act->y))!=NULL){//,x,y))!=NULL){
				//append the agent to nexthead
				A->AgentAppend(&(A->nexthead),c);
			}
			else{
				A->AgentFree(c);
				c=NULL;
			}
		}
		//TODO: check string lens of act and pass?


		//2: HEAL THE PARENT
		memset(act->f[act->ft],0,cpy*sizeof(char));

		csite->len = strlen(csite->S);

#ifdef DODEBUG
		if(csite->len==0){
			printf("zero length parent string!\n");
		}
#endif


		//Get rid of zero-length strings...
		//NB - grid status will be updated at the end of the timestep - simpler.
		if((dac = A->AgentCheckZeroLengthString(act))){
			//int x,y;
			switch(dac){
			case 1://Destroy active - only append passive
				A->AgentUnbind(pass,'P',1,act->spp,pass->spp);
				A->AgentAppend(&(A->nexthead),pass);
				//find_ag_gridpos(pass,run,&x,&y);
				//run->status[x][y]=G_NEXT;
				run->status[pass->x][pass->y]=G_NEXT;

				//find_ag_gridpos(act,run,&x,&y);
				run->grid[act->x][act->y]=NULL;
				run->status[act->x][act->y]=G_EMPTY;

				A->AgentFree(act);
				act = NULL;

				break;
			case 2://Destroy passive - only append active
				A->AgentUnbind(act,'A',1,act->spp,pass->spp);
				A->AgentAppend(&(A->nexthead),act);
				//find_ag_gridpos(act,run,&x,&y);
				run->status[act->x][act->y]=G_NEXT;


				//find_ag_gridpos(pass,run,&x,&y);
				run->grid[pass->x][pass->y]=NULL;
				run->status[pass->x][pass->y]=G_EMPTY;

				A->AgentFree(pass);
				pass = NULL;

				break;
			case 3://Destroy both
				printf("Destroying both parents after cleave!\nThis should never happen!\n");
				A->AgentUnbind(act,'A',1,act->spp,pass->spp);
				A->AgentUnbind(pass,'P',1,act->spp,pass->spp);
				A->AgentFree(act);
				act = NULL;
				A->AgentFree(pass);
				pass = NULL;
				break;
			default://This can't be right can it?
				if(act->ft == act->it){
					act->i[act->it]--;
				}
				break;
			}
		}
	}
	if(!dac){
		act->i[act->it]++;
	}

	return dac;
}





//todo(sjh): this should be part of the smspatial subclass
/*******************************************************************************
* @brief execute the current opcode in a reaction
*
* @param[in] A the stringmol bucket
*
* @param[in] run the grid information
*
* @param[in] act the first agent
*
* @param[in] pass the second agent
*
* @return 1 if decay happens, 0 if not
*******************************************************************************/
int ReactionExecuteOpcodeSpatial(stringPM *A, smsprun *run, s_ag *act, s_ag *pass){//, int x, int y){

	char *tmp;
	int safe_append=1;

	switch(*(act->i[act->it])){//*iptr[it]){

	case '$'://h-search
		//act->ft = act->it;
		char *cs;
		if(act->ft)
			cs = act->S;
		else
			cs = act->pass->S;
		tmp = OpcodeSearchInner(act->i[act->it],cs,A->blosum,&(act->it),&(act->ft),A->maxl);
		act->f[act->ft] = tmp;
		act->i[act->it]++;
		break;

	/*************
	 *   P_MOVE  *
	 *************/
	case '>':
			tmp=act->i[act->it];
			tmp++;
			switch(*tmp){
			case 'A':
				act->it = act->ft;
				act->i[act->it] = act->f[act->ft];
				act->i[act->it]++;
				break;
			case 'B':
				act->rt = act->ft;
				act->r[act->rt] = act->f[act->ft];
				act->i[act->it]++;
				break;
			case 'C':
				act->wt = act->ft;
				act->w[act->wt] = act->f[act->ft];
				act->i[act->it]++;
				break;
			default:
				act->it = act->ft;
				act->i[act->it] = act->f[act->ft];
				act->i[act->it]++;
				break;
			}
			break;


	/************
	 *   HCOPY  *
	 ************/
	case '='://h-copy
		//if(A->OpcodeCopy(act)<0){
		if(OpcodeCopy(act,A->domut,A->indelrate,A->subrate,A->maxl,
				A->blosum,A->granular_1,A->biomass)<0){
			A->AgentUnbind(act,'A',1,act->spp,pass->spp);
			A->AgentUnbind(pass,'P',1,act->spp,pass->spp);
		}
		break;


	/************
	 *   INC_R  *
	 ************/
	case '+'://h-copy
		if(A->granular_1==1){
			//printf("Incrementing read \n");
			/* Select the modifier */
			tmp=act->i[act->it];
			tmp++;
			switch(*tmp){
			case 'A':
				act->i[act->it]++;
				break;
			case 'B':
				act->r[act->rt]++;
				break;
			case 'C':
				act->w[act->wt]++;
				break;
			default:
				act->f[act->ft]++;
				break;
			}
		}
		act->i[act->it]++;
		break;



	/************
	 *  TOGGLE  *
	 ************/
	case '^'://p-toggle: toggle active pointer
			tmp=act->i[act->it];
			tmp++;
			switch(*tmp){
			case 'A':
				act->it = 1-act->it;
				break;
			case 'B':
				act->rt = 1-act->rt;
				break;
			case 'C':
				act->wt = 1-act->wt;
				break;
			default:
				act->ft = 1-act->ft;
				break;
			}
			act->i[act->it]++;
			break;

	/************
	 *  IFLABEL *
	 ************/
	case '?'://If-label
			act->i[act->it]=OpcodeIf(act->i[act->it],act->r[act->rt],act->S,A->blosum,A->maxl);
			break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			//Decide where to put the cleaved molecule
			if((/*dac = */OpcodeCleaveSpatial(A,run,act))){//,x,y))){
				//TODO: Need to determine what safe_append is used for (after looking at cleave)
				safe_append=0;	//extract_ag(&nowhead,p);
			}
			break;

	/**************
	 *  TERMINATE *
	 **************/
	case 0:
	case '}'://ex-end - finish execution

#ifdef V_VERBOSE
			printf("Unbinding...\n");
#endif

			A->AgentUnbind(act,'A',1,act->spp,pass->spp);
			run->status[act->x][act->y] = G_NEXT;

			A->AgentUnbind(pass,'P',1,act->spp,pass->spp);
			//int xx,yy;
			//find_ag_gridpos(pass,run,&xx,&yy);
			run->status[pass->x][pass->y] = G_NEXT;

			break;

	default://Just increment the i-pointer
		act->i[act->it]++;
		break;
	}
#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	print_exec(stdout,act,pass);
#endif


	//TODO: This action should be elsewhere - much harder to follow here
	if(safe_append){
		act->ect++;
		A->AgentAppend(&(A->nexthead),act);
		A->AgentAppend(&(A->nexthead),pass);
	}
	A->energy--;


	return 1;

}





//todo(sjh): this should be part of the smspatial subclass
/*******************************************************************************
* @brief attempt decay of a spatial agent
*
* @param[in] A the stringmol bucket
*
* @param[in] run the grid information
*
* @param[in] pag the agent
*
* @return 1 if decay happens, 0 if not
*******************************************************************************/
int AgentAttemptDecaySpatial(stringPM *A, smsprun *run, s_ag *pag){

 	float prob = A->decayrate;//1./pow(65,2);//4./3.); //This is now done in load_decay...

	float rno = RandomBetween0And1();

	if(rno<prob){
		//unbind_ag(pag);


		s_ag *bag;
		bag = NULL;//To prevend compiler "uninitialised" warning
		switch(pag->status){
		case B_UNBOUND:
			bag = NULL;
			break;
		case B_ACTIVE:
			bag = pag->pass;
			break;
		case B_PASSIVE:
			bag = pag->exec;
			break;
		}
		//int x,y;

		//find_ag_gridpos(pag,run,&x,&y);
		run->grid[pag->x][pag->y]=NULL;
		run->status[pag->x][pag->y]=G_EMPTY;

		A->AgentFree(pag);
		//TODO: sort this null-ing of free'd agents out!
		//pag = NULL;

		if(bag!=NULL){

			//find_ag_gridpos(bag,run,&x,&y);
			run->grid[bag->x][bag->y]=NULL;
			run->status[bag->x][bag->y]=G_EMPTY;

			A->AgentFree(bag);
			bag = NULL;
		}

		return 1;
	}
	else
		return 0;
}




/*******************************************************************************
* @brief get a filename for an input filename that doesn't exist
*
* @param[in] fn the filename of the file that might exist
*******************************************************************************/
void FilenameGetUnused(char *fn){

	int found = 1;
	char *point_pos, *tmp_pos;
	char tmp[128];
	char ofn[128];
	int ncopies=1;


	int ppos;
	strcpy(ofn,fn);

	point_pos = strchr(fn,'.');

	while(found){
		FILE * fpr;

		if((fpr=fopen(ofn,"r"))!=NULL){
			fclose(fpr);

			//Keep track of the original filename
			printf("Found a file called %s, incrementing counter...\n",ofn);
			/* STRATEGY: add a number at the file extension point
			 * until we find a file that hasn't been used yet...
			 */


			if(point_pos==NULL){
				printf("Can't extend the filename %s, exiting\n",fn);
				fflush(stdout);
				exit(37);
			}
			else{
				if(strlen(ofn)>127){
					printf("potential buffer overrun for filename %s, exiting\n",ofn);
					fflush(stdout);
					exit(38);
				}

				ppos = point_pos-fn;

				/*Build a new filename based on the original (ignoring the digits from intermediate attempts) */

				strncpy(tmp,fn,ppos);

				tmp_pos = &(tmp[ppos]);

				sprintf(tmp_pos,".%d%s",ncopies++,point_pos);

				strcpy(ofn,tmp);

			}
		}
		else{
			found = 0;
		}
	}
	//copy the result back to the filename
	strcpy(fn,ofn);
}





/*******************************************************************************
* @brief diagnostic printfs for MT load failure
*
* @details todo(sjh): "another attempt to do randseed properly"
*
* @param[in] fn the file name where the rng seed is found
*
* @param[in] printrandseed flag to print the seed
*
* @return the random seed
*******************************************************************************/
unsigned long RandomSeedInitFromFile(char *fn, int printrandseed=0){

	//TODO: more work needed with the rand seed and rng logic!!
	unsigned long seedin{42};
	unsigned int stmp{0};
	unsigned long sl{};

	bool foundrng{true};

	FILE *fpr;
	if((fpr=fopen(fn,"r"))!=NULL){
		//TODO: load the full RNG state using load_mt (RNGFILE in config)
		char *rngfn,*rngpath;
		rngfn=NULL;
		//rngfn = (char *) malloc (256*sizeof(char));
		rngpath = (char *) malloc (512*sizeof(char));
		memset(rngpath,0,512*sizeof(char));
		rngfn =  ParameterReadString(&fpr, "RNGFILE",0);

		if(rngfn !=NULL){
			//Add the file path from the input file to the RNGfile
			int l = strlen(fn);
			while(l>0){
				if(fn[l]=='/')
					break;
				else
					l--;
			}
			strncpy(rngpath,fn,l+1);
			char *pp;
			pp = &(rngpath[l+1]);
			strcpy(pp,rngfn);
			
			//TODO: check the logic of the following! Write in the usr guide..!
			FILE *fprng;
			if((fprng = fopen(rngpath,"r"))!=NULL){
				if(MersenneTwisterLoadState(rngpath) != load_mt_success){
					printf("ERROR reading Random Number Generator config %s\n",rngfn);
					foundrng = false;
				}
				else{
					//Record the RNG state (for debugging)
					FILE *rfp;
					rfp=fopen("RNGsmsp_initX_shouldwork.dat","w");
                    MersenneTwisterPrintStatusToFile(rfp);
					fclose(rfp);
				}
				fclose(fprng);
			}
			else{
				printf("ERROR reading Random Number Generator config %s\n",rngfn);
				foundrng = false;
			}
		}
		else{
			foundrng=false;
		}

		int rerr = ParameterReadOrDefineUnsignedInt(fn,"RANDSEED", &stmp, seedin, 0);//read_param_int(fpr,"RANDSEED",&stmp,1);


		if(rerr){
			printf("Error %d reading RANDSEED\n",rerr);
			exit(0);
		}

		free(rngfn);
		free(rngpath);
		fclose(fpr);
	}


	unsigned long rseed;
	if(!foundrng){
		if(stmp){//This means we have read it from the file...
			sl = stmp;
			rseed = RandomInitLong(&sl);
		}
		else{
			rseed = RandomInitLong(NULL);
		}
	}
	else{
		rseed = seedin = sl = stmp;
	}

	if(printrandseed){
		char frfn[128];

		sprintf(frfn,"randseed.txt");

		FilenameGetUnused(&(frfn[0]));

		FILE *frs;
		if((frs=fopen(frfn,"w"))==NULL){
			printf("Coundln't open %s\n",frfn);
			getchar();
			exit(39);
		}else{
			fprintf(frs,"(unsigned) random seed is %lu (%lu was seedin)  \n",rseed,seedin);
			fflush(frs);
			fclose(frs);
		}
	}

	return rseed;
}





//todo(sjh): integrate this with stringPM::TimestepUpdate
/*******************************************************************************
* @brief diagnostic printfs for MT load failure
*
* @details see Mersenne Twister documentation
*
* @param[in] ec error code
*
* @param[in] posn position of the error in the file
*
* @param[in] fn the file name
*
* @return mt_error_code
*******************************************************************************/
int TimestepIncrementSpatial(stringPM *A, smsprun *run){

	//again, we follow TimestepIncrement, but are a little more careful with the binding and uncoupling
	s_ag *pag;

	//A->energy += A->estep;

	//Set the energy to the max possible number of molecules...
	//printf("Before update, energy is %d\n",A->energy);
	A->energy = run->gridx * run->gridy;
	//printf("After update, energy is %d\n",A->energy);

	int ct=0;

	while(A->nowhead!=NULL){
 
 		s_ag *bag;
		pag = A->AgentSelectRandomly(A->nowhead,-1);
		A->AgentExtract(&A->nowhead,pag);

		//For debugging RNG diffs.
		if(A->timestep == 90001){
			printf("%d\t%d\t%d\t\n",//%c%c%c%c\n",
					ct++,
					RandomNumberGeneratorGetState(),
					pag->idx//,
					//&pag->i[pag->it],
					//&pag->r[pag->rt],
					//&pag->w[pag->wt],
					//&pag->f[pag->ft]
					);
		}

		//extract any partner:
		bag = NULL;

		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			bag = pag->pass;
			A->AgentExtract(&(A->nowhead),bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			A->AgentExtract(&(A->nowhead),bag);
			break;
		}

		if(!AgentAttemptDecaySpatial(A,run,pag)){
			int changed = 0;
			if(A->energy>0){
				switch(pag->status){
				case B_UNBOUND:
					//seek binding partner, set binding states.
					//changed = A->testbind(pag);

					//int x,y;
					align sw;

					//We can only bind neighbours in the spatial model
					//find_ag_gridpos(pag,run,&x,&y);
					run->status[pag->x][pag->y]=G_NEXT;

					if((bag = ReactionSeekRandomSpatialPartner(A,run,pag->x,pag->y))!=NULL){

						//TODO: We need to make sure that bag is in nowhead first!
						A->AgentExtract(&(A->nowhead),bag);

						//Now we've found a potential partner, we can see if it binds:
						float bprob;
						bprob = A->AgentsAlign(pag,bag,&sw);

						float rno;
						rno = RandomBetween0And1();
						if(rno<bprob){//Binding success!
							//figure out which is the executing string:
							A->ReactionSetupExecution(pag,bag,&sw);
							pag->nbind++;
							bag->nbind++;

							A->energy--;

							A->AgentAppend(&(A->nexthead),pag);
							A->AgentAppend(&(A->nexthead),bag);
							changed=1;
						}
					}

					break;
				case B_PASSIVE:


					//find_ag_gridpos(pag->exec,run,&x,&y);

					changed = ReactionExecuteOpcodeSpatial(A,run,pag->exec,pag);//,x,y);

					break;
				case B_ACTIVE:

					//find_ag_gridpos(pag,run,&x,&y);
					changed = ReactionExecuteOpcodeSpatial(A,run,pag,pag->pass);//,x,y);
					break;
				default:
					printf("ERROR: agent with unknown state encountered!\n");
				}
			}
			if(!changed){
				A->AgentAppend(&(A->nexthead),pag);
				if(bag!=NULL)
					A->AgentAppend(&(A->nexthead),bag);

			}
		}
	}

	A->UpdateNowNext();
	TimestepGridIncrement(run);



	//#ifdef DEBUG
	for(int x=0;x<run->gridx;x++){
		for(int y=0;y<run->gridy;y++){
			if(run->grid[x][y]!=NULL){
				if(!(A->AgentAddressInList(A->nowhead,run->grid[x][y]))){
					printf("Agent species %d not in nowhead at %d, %d\n", run->grid[x][y]->spp->spp, x, y);
				}
			}
		}
	}

	return 0;
}





//todo(sjh): integrate the smsprun struct with stringPM
/*******************************************************************************
* @brief reads a config file and sets up a Spatial Stringmol run
*
* @param[in] fn the config file name
*
* @param[in] A the stringPM object
*
* @param[in] run grid data in an smsprun struct
*
* @param[in] runno the run number
*
* @return 0 always
*******************************************************************************/
int StringmolSpatialConfigureFromFile(const char *fn, stringPM *A, smsprun **run, int runno){

	A->ParametersLoad(fn,0,1);
	A->AgentsLoad(fn,NULL,0,0);

	A->run_number = runno;

	//Initialize the popdy file...
	PopdyInitFile(A);

	//TODO(sjh): It's a little perverse getting this run object out, but we have to decide whether the grid is 'core' stringmol...
	*run = A->grid;
	//if(run == NULL){
	if(*run == NULL){
		printf("No grid data entered for spatial stringmol!\nexiting...");
		exit(33);
	}


	//TODO: Now we have to place each agent on the grid - use the makenext() model -
	//But we only need to do this if extit == 0, because otherwise we'll have the molecular positions...
	if(!A->timestep){
		while(A->nowhead!=NULL){
			s_ag *pag;
			pag = A->AgentSelectRandomly(A->nowhead,-1);
			A->AgentExtract(&(A->nowhead),pag);
			int found = 0;

			//Check to see if a position has been set for each molecule..
			if( pag->x > -1 ){
				if( pag->y > -1){
					AgentPlaceOnGrid(pag,*run,pag->x,pag->y);
					A->AgentAppend(&(A->nexthead),pag);
					found = 1;
				}
				else{
					printf("Bad xy position for this agent (%d,%d) \n",pag->x,pag->y);
					exit(0);
				}

			}

			while(!found){
				int pos = (*run)->gridx * (*run)->gridy * RandomBetween0And1();

				int x = pos%(*run)->gridx;
				int y = pos/(*run)->gridx;

				if((*run)->grid[x][y]==0){

					//TODO: this command should be moved to the smspatial
					//((uint8_t *)screen->pixels)[x + (y * sdlPitch)] = 0;

					//Add the partner in the Moore neighborhood
					int ffound=0;

					while(!ffound){
						int xx,yy;
						//randy_Moore(const int X, const int Y, const int Xlim, const int Ylim, int *xout, int *yout){
						GridSelectRandomMooreNeighbour(x,y,(*run)->gridx,(*run)->gridy,&xx,&yy);
						if((*run)->grid[xx][yy]==0){

							ffound=found=1;

							//Place each agent on the list
							AgentPlaceOnGrid(pag,*run,x,y);

							//Move to the 'used' bucket
							A->AgentAppend(&(A->nexthead),pag);

							s_ag *bag;
							bag = A->AgentSelectRandomly(A->nowhead,-1);
							if(bag != NULL){
								A->AgentExtract(&(A->nowhead),bag);
								AgentPlaceOnGrid(bag,*run,xx,yy);
								A->AgentAppend(&(A->nexthead),bag);
							}
						}
					}
				}
			}
		}

		A->UpdateNowNext();
	}
	else{
		//Fill the grid
		s_ag *pag{};

		for(int x=0;x<(*run)->gridx;x++)
			for(int y=0;y<(*run)->gridy;y++)
				(*run)->grid[x][y] =NULL;

		for(pag=A->nowhead; pag != NULL; pag=pag->next){
			AgentPlaceOnGrid(pag,*run,pag->x,pag->y);
		}

		TimestepGridIncrement(A->grid);
	}


	//#ifdef DEBUG
	for(int x=0;x<(*run)->gridx;x++){
		for(int y=0;y<(*run)->gridy;y++){
			if((*run)->grid[x][y]!=NULL){
				if(!(A->AgentAddressInList(A->nowhead,(*run)->grid[x][y]))){
					printf("Agent species %d not in nowhead at %d, %d\n", (*run)->grid[x][y]->spp->spp, x, y);
				}
			}
		}
	}
	//#endif

	return 0;
}





/*******************************************************************************
* @brief stringmol on a grid
*
* @details "Nothing makes sense in evolutin except in the light of parasitism"
*
* @param[in] argc number of arguments
*
* @param[in] argv the arguments
*
* @return 0 unless there's an error
*******************************************************************************/
int StringmolSpatial(int argc, char *argv[]) {

	printf("Hello spatial stringmol world\n");

	SMspp		SP{};
	stringPM	A{&SP};

	smsprun *run{};
	run = NULL;

	A.randseed = RandomSeedInitFromFile(argv[2]);
	StringmolSpatialConfigureFromFile(argv[2],&A,&run,1);

	int bt,ct{0};
	ct = A.AgentsCount(A.nowhead,-1);
	printf("Initialisation done, number of molecules is %d\n",ct);

	//This used to be called here - but better to do it before smspatial_init()
	//that way you don't reinitialize the seed after (perhaps) positioning the
	//agents on the grid
	//
	//A.randseed = init_randseed(argv[2]);


	//reload file: should have identical settings to the input conig:
	FILE *fpp{};
	char fn[128]{};
	sprintf(fn,"reload_%05u.conf",A.timestep);
	fpp = fopen(fn,"w");
	A.print_conf(fpp);
	fclose(fpp);


	//Graphic to console - sometimes useful..
	//A.print_grid(stdout);


//	while(A.nagents(A.nowhead,-1)){
	while((A.timestep < A.nsteps) && (ct > 0)){

		//if(!(A.extit%100) || A.extit==1){
		//if(!(A.extit%100)){
		if(!(A.timestep%A.image_every)){
			bt = ct - A.AgentsCount(A.nowhead,B_UNBOUND);
			printf("Step %u done, number of molecules is %d, nbound = %d\n",A.timestep,ct,bt);

			GridSavePNG(&A, smpic_spp);
			GridSavePNG(&A, smpic_len);
		}

		//When debugging you can set specific timesteps here (not elegant, but quick!)
		//if(!(A.extit%1000) || A.extit==1
		//		|| (A.extit%1000) == 999 || (A.extit%1000) == 1  || !(A.extit%97) || !(A.extit%47)
		//		|| (A.extit%1000) == 99 || (A.extit%1000) == 100|| (A.extit%1000) == 101
		//		|| (A.extit>90001 && A.extit <9000))
		if(!(A.timestep%A.report_every)){
			
			
			sprintf(fn,"out1_%05u.conf",A.timestep);
			fpp = fopen(fn,"w");

			A.print_conf(fpp);
			fclose(fpp);

			FILE *fp{};
			sprintf(fn,"splist%u.dat",A.timestep);
			fp = fopen(fn,"w");
			SP.SpeciesListPrint(fp);
			fclose(fp);


			SpeciesPrintCounts(&A,A.timestep);

		}

		TimestepIncrementSpatial(&A,run);

#ifdef DODEBUG
		printf("Nowhead is %p, Nexthead is %p\n",A.nowhead,A.nexthead);
		s_ag *p;
		p=A.nowhead;
		int mno=0;
		while(p!=NULL){
			int x,y;
			find_ag_gridpos(p,run,&x,&y);
			printf("%d, %p, [%d,%d]  status: %d, bound to %p / %p, prev = %p, next = %p\n",
			++mno,p,x,y,p->status,p->exec,p->pass,p->prev,p->next);
			p = p->next;
		}
#endif

		A.timestep++;
		ct = A.AgentsCount(A.nowhead,-1);

	}

	printf("FINISHED smspatial\n");
	fflush(stdout);
	return 0;
}





/*******************************************************************************
* @brief wrapper function for lodepng to write png to file
*
* @details Encode from raw pixels to disk with a single function call. The image
*          argument has width * height RGBA pixels or width * height * 4 bytes
*
* @param[in] filename the png file name
*
* @param[in] image the image data
*
* @param[in] width the width of the image (x)
*
* @param[in] height the height of the image (y)
*
* @return 0 always
*******************************************************************************/
void PNGEncodeAndSave(const char* filename, const std::vector<unsigned char>& image, unsigned width, unsigned height)
{
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}





/*******************************************************************************
* @brief save the grid state as a png
*
* @details "Nothing makes sense in evolutin except in the light of parasitism"
*
* @param[in] A the stringmol bucket
*
* @param[in] pt enum to say whether length or species no. will be coloured.
*
* @return 0 always
*******************************************************************************/
int GridSavePNG(stringPM *A, smpic pt){


	smsprun *run;
	int colours[70][3] = {
	{ 0, 0, 167 },	{ 0, 10, 177 },	{ 0, 20, 186 },	{ 0, 29, 196 },	{ 0, 39, 206 },	{ 0, 49, 216 },	{ 0, 59, 226 },	{ 0, 69, 235 },
	{ 0, 78, 245 },	{ 0, 88, 255 },	{ 0, 98, 255 },	{ 0, 108, 255 },{ 0, 118, 255 },{ 0, 128, 255 },{ 0, 137, 255 },{ 0, 147, 255 },
	{ 0, 157, 255 },{ 0, 167, 255 },{ 10, 177, 255 },{ 20, 186, 255 },{ 29, 196, 255 },{ 39, 206, 255 },{ 49, 216, 255 },{ 59, 226, 255 },
	{ 69, 235, 255 },{ 78, 245, 255 },{ 88, 255, 255 },	{ 98, 255, 245 },{ 108, 255, 235 },	{ 118, 255, 226 },	{ 128, 255, 216 },	{ 137, 255, 206 },
	{ 147, 255, 196 },{ 157, 255, 186 },{ 167, 255, 177 },{ 177, 255, 167 },{ 186, 255, 157 },{ 196, 255, 147 },{ 206, 255, 137 },	{ 216, 255, 128 },{ 226, 255, 118 },
	{ 235, 255, 108 },	{ 245, 255, 98 },{ 255, 255, 88 },{ 255, 245, 78 },	{ 255, 235, 69 },	{ 255, 226, 59 },	{ 255, 216, 49 },	{ 255, 206, 39 },
	{ 255, 196, 29 },	{ 255, 186, 20 },	{ 255, 177, 10 },	{ 255, 167, 0 },	{ 255, 157, 0 },	{ 255, 147, 0 },	{ 255, 137, 0 },
	{ 255, 128, 0 },	{ 255, 118, 0 },	{ 255, 108, 0 },	{ 255, 98, 0 },	{ 255, 88, 0 },	{ 245, 78, 0 },	{ 235, 69, 0 },	{ 226, 59, 0 },
	{ 216, 49, 0 },	{ 206, 39, 0 },	{ 196, 29, 0 },	{ 186, 20, 0 },	{ 177, 10, 0 },	{ 167, 0, 0 },
	};




	run = A->grid;

	//Create the PNG
	std::vector<unsigned char> image;
	image.resize(run->gridx * run->gridy * 4);

	int x,y,val;

	//Bit masks for the 8*8*4 rgb cube
	int rmask = 0b11100000;
	int gmask = 0b00011100;
	int bmask = 0b00000011;
	int len;


	for(x=0;x<run->gridx;++x){
		for (y=0;y<run->gridy;++y) {

			if(run->grid[x][y] == NULL){

				image[4 * run->gridx * y + 4 * x + 0] = 0;//255 * !(x & y);
				image[4 * run->gridx * y + 4 * x + 1] = 0;//x ^ y;
				image[4 * run->gridx * y + 4 * x + 2] = 0;//x | y;
				image[4 * run->gridx * y + 4 * x + 3] = 255;



			}
			else{
				switch(pt){
				case smpic_len:
					len = strlen(run->grid[x][y]->spp->S);

					if(len>69){
						image[4 * run->gridx * y + 4 * x + 0] = 255;//255 * !(x & y);
						image[4 * run->gridx * y + 4 * x + 1] = 255;//x ^ y;
						image[4 * run->gridx * y + 4 * x + 2] = 255;//x | y;
						image[4 * run->gridx * y + 4 * x + 3] = 255;

					}
					else{

						image[4 * run->gridx * y + 4 * x + 0] = colours[len][0];//255 * !(x & y);
						image[4 * run->gridx * y + 4 * x + 1] = colours[len][1];//x ^ y;
						image[4 * run->gridx * y + 4 * x + 2] = colours[len][2];//x | y;
						image[4 * run->gridx * y + 4 * x + 3] = 255;

					}
					break;
				case smpic_spp:

					val = ((run->grid[x][y]->spp->spp)*5) % 256;

					image[4 * run->gridx * y + 4 * x + 0] = (32  * (1+((val & rmask) >> 5)))-1;
					image[4 * run->gridx * y + 4 * x + 1] = (32  * (1+((val & gmask) >> 2)))-1;
					image[4 * run->gridx * y + 4 * x + 2] = (64 *  (1+((val & bmask)     )))-1;
					image[4 * run->gridx * y + 4 * x + 3] = 255;

					break;
				}

			}
		}
	}


	char filename[128];
	switch(pt){
	case smpic_len:
		sprintf(filename,"lenframe%07u.png",A->timestep);
		break;
	case smpic_spp:
		sprintf(filename,"sppframe%07u.png",A->timestep);
		break;

	}
	PNGEncodeAndSave(filename, image, run->gridx, run->gridy);

	return 0;
}





//todo(sjh): The fails with segfaults - fix!
/*******************************************************************************
* @brief stringmol on a grid - calculate ancestry
*
* @details Strategy is to create a text file containing the ancestry. We'll
*          write an R or graphviz script to parse this and generate figures.
*          argv[1] 35
*          argv[2] first timestep for pics
*          argv[3] last timestep for pics
*          argv[4] step between timesteps
*
* @param[in] argc number of arguments
*
* @param[in] argv the arguments
*
* @return 0 always
*******************************************************************************/
int StringmolSpatialPicsFromLogs(int argc, char *argv[]){

	/*This colormap was generated using the following R script:
	 *
			require(colorRamps)
			x <- matlab.like(70)
			y <- col2rgb(x)
			dim <- 70
			message(sprintf("int colours[%d][3] = {",dim))
			for(i in 0:dim) message(sprintf("{ %d, %d, %d },",y[1,i],y[2,i],y[3,i]))
			message("};")
	 *
	 */


	int colours[70][3] = {
	{ 0, 0, 167 },	{ 0, 10, 177 },	{ 0, 20, 186 },	{ 0, 29, 196 },	{ 0, 39, 206 },	{ 0, 49, 216 },	{ 0, 59, 226 },	{ 0, 69, 235 },
	{ 0, 78, 245 },	{ 0, 88, 255 },	{ 0, 98, 255 },	{ 0, 108, 255 },{ 0, 118, 255 },{ 0, 128, 255 },{ 0, 137, 255 },{ 0, 147, 255 },
	{ 0, 157, 255 },{ 0, 167, 255 },{ 10, 177, 255 },{ 20, 186, 255 },{ 29, 196, 255 },{ 39, 206, 255 },{ 49, 216, 255 },{ 59, 226, 255 },
	{ 69, 235, 255 },{ 78, 245, 255 },{ 88, 255, 255 },	{ 98, 255, 245 },{ 108, 255, 235 },	{ 118, 255, 226 },	{ 128, 255, 216 },	{ 137, 255, 206 },
	{ 147, 255, 196 },{ 157, 255, 186 },{ 167, 255, 177 },{ 177, 255, 167 },{ 186, 255, 157 },{ 196, 255, 147 },{ 206, 255, 137 },	{ 216, 255, 128 },{ 226, 255, 118 },
	{ 235, 255, 108 },	{ 245, 255, 98 },{ 255, 255, 88 },{ 255, 245, 78 },	{ 255, 235, 69 },	{ 255, 226, 59 },	{ 255, 216, 49 },	{ 255, 206, 39 },
	{ 255, 196, 29 },	{ 255, 186, 20 },	{ 255, 177, 10 },	{ 255, 167, 0 },	{ 255, 157, 0 },	{ 255, 147, 0 },	{ 255, 137, 0 },
	{ 255, 128, 0 },	{ 255, 118, 0 },	{ 255, 108, 0 },	{ 255, 98, 0 },	{ 255, 88, 0 },	{ 245, 78, 0 },	{ 235, 69, 0 },	{ 226, 59, 0 },
	{ 216, 49, 0 },	{ 206, 39, 0 },	{ 196, 29, 0 },	{ 186, 20, 0 },	{ 177, 10, 0 },	{ 167, 0, 0 },
	};

	////generates pics from a config file via SDL
	//void pics_from_config(int s, int f, int step){

	int s = atoi(argv[2]);
	int f = atoi(argv[3]);
	int step = atoi(argv[4]);


	char fn[256];

	SMspp		SP;
	stringPM	A(&SP);

	smsprun *run;

	StringmolSpatialConfigureFromFile("out1_00100.conf",&A,&run,1);


	for(int i=s;i<=f;i+=step){

		memset(fn,0,256*sizeof(char));

		sprintf(fn,"out1_%05d.conf",i);


		if((StringmolSpatialConfigureFromFile(fn,&A,&run,1))==0){

			//Should be able to dump the grid image now...

			//Create the PNG
			std::vector<unsigned char> image;
			image.resize(run->gridx * run->gridy * 4);

			int x,y;

			if(i==s){
				for(x=0;x<run->gridx;++x){
					for (y=0;y<run->gridy;++y) {

						if(y>69){
							image[4 * run->gridx * y + 4 * x + 0] = 255;//255 * !(x & y);
							image[4 * run->gridx * y + 4 * x + 1] = 255;//x ^ y;
							image[4 * run->gridx * y + 4 * x + 2] = 255;//x | y;
							image[4 * run->gridx * y + 4 * x + 3] = 255;

						}
						else{

							image[4 * run->gridx * y + 4 * x + 0] = colours[y][0];//255 * !(x & y);
							image[4 * run->gridx * y + 4 * x + 1] = colours[y][1];//x ^ y;
							image[4 * run->gridx * y + 4 * x + 2] = colours[y][2];//x | y;
							image[4 * run->gridx * y + 4 * x + 3] = 255;

						}

					}
				}

				char filename[128];
				sprintf(filename,"lenkey.png");
				PNGEncodeAndSave(filename, image, run->gridx, run->gridy);


				/* Make the rgb species pallette
				 *
				 */
				//Bit masks for the 8*8*4 rgb cube
				int rmask = 0b11100000;
				int gmask = 0b00011100;
				int bmask = 0b00000011;

				
				printf("Mask values are r:%d, g:%d, b:%d\n",rmask,gmask,bmask);

				for(y=0;y<run->gridy;++y){
					int val = y;//((x+1) * 5) % 256;

					int rbits = (val & rmask);
					int gbits = (val & gmask);
					int bbits = (val & bmask);

					printf("y = %d, val = %d, \tRGB unscaled = r:%d, g:%d, b:%d \t scaled = r:%d g:%d b:%d\n",y,val,rbits,gbits,bbits,rbits>>5,gbits>>2,bbits);


				}
				for(x=0;x<run->gridx;++x){
					for (y=0;y<run->gridy;++y) {

						int val = ((y+1)*5) % 256;

						image[4 * run->gridx * y + 4 * x + 0] = (32  * (1+((val & rmask) >> 5)))-1;
						image[4 * run->gridx * y + 4 * x + 1] = (32  * (1+((val & gmask) >> 2)))-1;
						image[4 * run->gridx * y + 4 * x + 2] = (64 *  (1+((val & bmask)     )))-1;
						image[4 * run->gridx * y + 4 * x + 3] = 255;


					}
				}

				sprintf(filename,"sppkey.png");
				PNGEncodeAndSave(filename, image, run->gridx, run->gridy);

			}


			for(x=0;x<run->gridx;++x){
				for (y=0;y<run->gridy;++y) {


					if(run->grid[x][y] == NULL){

						image[4 * run->gridx * y + 4 * x + 0] = 0;//255 * !(x & y);
						image[4 * run->gridx * y + 4 * x + 1] = 0;//x ^ y;
						image[4 * run->gridx * y + 4 * x + 2] = 0;//x | y;
						image[4 * run->gridx * y + 4 * x + 3] = 255;



					}
					else{

						int len = strlen(run->grid[x][y]->spp->S);

						if(len>69){
							image[4 * run->gridx * y + 4 * x + 0] = 255;//255 * !(x & y);
							image[4 * run->gridx * y + 4 * x + 1] = 255;//x ^ y;
							image[4 * run->gridx * y + 4 * x + 2] = 255;//x | y;
							image[4 * run->gridx * y + 4 * x + 3] = 255;

						}
						else{

							image[4 * run->gridx * y + 4 * x + 0] = colours[len][0];//255 * !(x & y);
							image[4 * run->gridx * y + 4 * x + 1] = colours[len][1];//x ^ y;
							image[4 * run->gridx * y + 4 * x + 2] = colours[len][2];//x | y;
							image[4 * run->gridx * y + 4 * x + 3] = 255;

						}
					}
				}
			}


			char filename[128];
			sprintf(filename,"lenframe%07u.png",A.timestep);
			PNGEncodeAndSave(filename, image, run->gridx, run->gridy);

			A.BucketReset();
			SP.SpeciesListClear();


		}
		else{
			printf("File %s not found, exiting\n",fn);
			exit(34);
		}
	}


	return 0;
}





struct anc_node{
	int spno;
	int origintime;
	char * S;
	anc_node *aa; //active parent
	anc_node *pp; //passive parent
	anc_node *dd; //pointer to descendents (may be useful...)
	anc_node *ss; //pointer to siblings (may be useful...)
};





/*******************************************************************************
* @brief allocate memory for a new ancestry node
*
* @return the node
*******************************************************************************/
anc_node * AncestryAllocNode(){
	anc_node *node;

	node = static_cast<anc_node *>(malloc(sizeof(anc_node)));

	node->spno=0;
	node->origintime=0;

	node->S = NULL;

	node->aa = NULL;
	node->pp = NULL;
	node->dd = NULL;
	node->ss = NULL;

	return node;
}





/*******************************************************************************
* @brief create a new ancestry node
*
* @param[in] sp the species
*
* @param[in] time the time found
*
* @return 0 always
*******************************************************************************/
anc_node * AncestryNewNode(int sp, int time){
	anc_node *node;

	//node = AncestryAllocNode();
	node = static_cast<anc_node *>(malloc(sizeof(anc_node)));

	node->spno=0;
	node->origintime=0;

	node->S = NULL;

	node->aa = NULL;
	node->pp = NULL;
	node->dd = NULL;
	node->ss = NULL;


	node->spno = sp;
	node->origintime = time;
	return node;
}





/*******************************************************************************
* @brief find the parents of a species
*
* @param[in] aa anc_node - position in the ancestry
*
* @param[in] timestep
*
* @param[in] depth how far back to search
*
* @param[in] ofn output file name
*
* @param[in] found_spp the number of species found so far
*
* @return 0 always
*******************************************************************************/
void AncestryFindParents(anc_node *aa, int timestep, int depth, char * ofn, int *found_spp){

	FILE *ofp;
	const int llen = 3000;

	char fn[128],line[llen],seq[llen];

	int s1,a1,p1;
	int nr,t1;

	int found = 0;

	while(!found){
		sprintf(fn,"splist%d.dat",timestep);
		FILE *fp;
		if((fp=fopen(fn,"r"))==NULL){
			printf("Cannot open file %s so cannot proceed\n",fn);
			return;
		}
		else{

			timestep -= 100;
			int lno = 0;
			printf("Successfully opened file %s\nNow looking for species %d\n", fn, aa->spno);
			while((fgets(line,llen,fp)) != NULL){

				//56189,-1,-1,1,EK$OYHOJLRWEK$BLUBO^>C$=?>D$BLUB}B$=?$$$BLBLUB}OYYHKOBLBLUB}OYYHKO
				sscanf(line,"%d,%d,%d",&s1,&a1,&p1);
				lno++;

				//TODO: Need a better way to resolve multiple parenting events, especially when spread over many logfiles....
				if(s1 == aa->spno && s1>a1 && s1 >p1){
					if(a1==-1 && p1==-1){
						//printf("Line %d: Ancestors of spp %d are from an earlier file, decrementing timestep\n",lno,s1);
						//fclose(fp);
						if(!timestep){
							printf("Line %d: Made it to the dawn of time!\nLUCA is:\n",lno);
							//1,-1,-1,1,WWGEWLHHHRLUEUWJJJRJXUUUDYGRHJLRWWRE$BLUBO^B>C$=?>$$BLUBO%}OYHOB
							sscanf(line,"%*d,%*d,%*d,%*d,%2000s",&(seq[0]));
							printf("depth\ttime\tspp\tact\tpass\tsequence\n");
							printf("%d\t1\t%d\t%d\t%d\t%s\n",depth,s1,a1,p1,seq);

							//Add this discovery to the ancestry log file
							ofp=fopen(ofn,"a");

							//fprintf(ofp,"depth,time,spp,act,pass,sequence\n");
							//todo: figure out how graphviz will view this information
							fprintf(ofp,"%d,1,%d,-1,-1,%s\n",depth,aa->spno,seq);
							fclose(ofp);

							fclose(fp);
							return;
						}
						break;
					}
					else{
						printf("Line %d: Ancestors of spp %d found! Active is %d, Passive is %d\n",lno,s1,a1,p1);
						sscanf(line,"%*d,%*d,%*d,%*f,%d,%d,%*d,%2000s",&nr,&t1,&(seq[0]));
						found = 1;
						//fclose(fp);


						//Add this discovery to the ancestry log file
						ofp=fopen(ofn,"a");

						//fprintf(ofp,"depth,time,spp,act,pass,sequence\n");
						fprintf(ofp,"%d,%d,%d,%d,%d,%s\n",depth,t1,aa->spno,a1,p1,seq);

						printf("depth\ttime\tspp\tact\tpass\tsequence\n");
						printf("%d\t%d\t%d\t%d\t%d\t%s\n",depth,t1,s1,a1,p1,seq);


						fclose(ofp);
						break;
					}
				}
			}
			//if(fp != NULL)
			//fclose(fp);		printf("Finished looking at file %s\n");
			fclose(fp);

		}
	}


	aa->origintime = t1;

	int L = strlen(seq)+1;
	aa->S = (char *) malloc( L *sizeof(char));
	memset(aa->S,0,L*sizeof(char));
	strcpy(aa->S,seq);

	//Go up a generation:

	if(!found_spp[aa->spno]){

		found_spp[aa->spno]=1;
		aa->aa = AncestryNewNode(a1,-1);
		AncestryFindParents(aa->aa,timestep,depth+1,ofn,found_spp);
		//Don't bother tracing the passive branch if parents are the same species..
		if(a1 != p1){
			aa->pp = AncestryNewNode(p1,-1);
			AncestryFindParents(aa->pp,timestep,depth+1,ofn,found_spp);
		}
	}

}






/*******************************************************************************
* @brief stringmol on a grid - calculate ancestry
*
* @details Strategy is to create a text file containing the ancestry. We'll
*          write an R or graphviz script to parse this and generate figures.
*          argv[1] 34
*          argv[2] species number
*          argv[3] timestep
*          argv[4] outfile name
*
* @param[in] argc number of arguments
*
* @param[in] argv the arguments
*
* @return 0 always
*******************************************************************************/
int StringmolSpatialAncestry(int argc, char *argv[]){


	if(argc != 4 && argc !=5){
		printf("Try again :) Usage:\n\tstringmol 34 spno timestep (outfile)\n");
		exit(340);
	}

	char ofn[128];
	FILE *ofp;

	if(argc == 5){
		sprintf(ofn,"%s",argv[4]);
	}
	else{
		sprintf(ofn,"outfile.txt");
	}

	if((ofp = fopen(ofn,"w"))==NULL){
		printf("Problem opening output file %s\n",ofn);
		exit(341);
	}

	fprintf(ofp,"depth,time,spp,act,pass,sequence\n");
	fclose(ofp);

	anc_node *aa;

	aa = AncestryAllocNode();

	int spno = atoi(argv[2]);
	int timestep = atoi(argv[3]);
	//const int llen = 3000;

	int *found_spp;

	found_spp = (int *)malloc((spno+1)*sizeof(int));
	memset(found_spp,0,(spno+1)*sizeof(int));

	/*set up the recursion*/

	aa->spno = spno;

	AncestryFindParents(aa,timestep,0,ofn,found_spp);

	/*
	while(!found_luca){

		sprintf(fn,"splist%d.dat",timestep);
		if((fp=fopen(fn,"r"))==NULL){
			printf("Cannot open file %s so cannot proceed\n",fn);
			exit(341);
		}
		else{

			printf("Successfully opened file %s\nNow looking for species %d\n", fn, spno);
			while((fgets(line,llen,fp)) != NULL){

				//56189,-1,-1,1,EK$OYHOJLRWEK$BLUBO^>C$=?>D$BLUB}B$=?$$$BLBLUB}OYYHKOBLBLUB}OYYHKO
				sscanf(line,"%d,%d,%d",&s1,&a1,&p1);

				if(s1 == spno){
					if(a1==-1 && p1==-1){
						printf("Ancestors are from earlier run, decrementing timestep\n");
						timestep -= 100;
						fclose(fp);
						break;
					}
					else{
						printf("Ancestors found! Active is %d, Passive is %d\n",a1,p1);
						break;
					}
				}
			}

			//find the species
		}
	}
	*/

	/*Open the splist file and find the entries for the species*/

	return 0;
}





struct comm_node{
	int aspno;
	int pspno;
	int count;
	char * Sa;
	char * Sp;
	comm_node *next;
};





/*******************************************************************************
* @brief Allocate and initialise a new comm_node object
*
* @return the node object. todo(sjh): no failure reporting
*******************************************************************************/
comm_node * ReactionObservedAlloc(){
	comm_node *node;

	node = static_cast<comm_node *>(malloc(sizeof(comm_node)) );

	node->aspno=0;
	node->pspno=0;
	node->count=0;

	node->Sa = NULL;
	node->Sp = NULL;

	node->next = NULL;
	return node;
}





/*******************************************************************************
* @brief allocate memory and copy the string into it
*
* @param[in] S the sequence string
*
* @return the allocated string
*******************************************************************************/
char * SequenceAllocAndCopy(char *S){
	char *out;
	int len = 1+strlen(S);
	out = (char *) malloc((len)*sizeof(char));
	memset(out,0,len*sizeof(char));
	//strncpy(out,S,(len-1));
	strcpy(out,S);
	return out;

}





/*******************************************************************************
* @brief A new comm_node object recording an observed bind
*
* @param[in] ac the active species number
*
* @param[in] pa the passive species number
*
* @param[in] Sa the active species' sequence
*
* @param[in] Sp the passive species' sequence
*
* @return the node
*******************************************************************************/
comm_node * ReactionNewObservedPairing(int ac, int pa, char *Sa, char *Sp){
	comm_node *aa;

	aa = ReactionObservedAlloc();

	aa->aspno = ac;
	aa->pspno = pa;
	aa->count = 1;

	aa->Sa = SequenceAllocAndCopy(Sa);
	aa->Sp = SequenceAllocAndCopy(Sp);


	//TODO: copy the strings - if necessary...

	return aa;
}





//todo(sjh) have a 1m timestep file to see what this puts out!
/*******************************************************************************
* @brief stringmol on a grid - calculate community - not sure how this works!
*
* @details Strategy is to create a text file containing the ancestry. We'll
*          write an R or graphviz script to parse this and generate figures.
*          argv[1] 36
*          argv[2] config file
*
* @param[in] argc number of arguments
*
* @param[in] argv the arguments
*
* @return 0 always
*******************************************************************************/
int StringmolSpatialCommunity(int argc, char *argv[]){

	comm_node * bbb;
	bbb = NULL;

	if(argc != 3){
		printf("Try again :) Usage:\n\tstringmol 36 outfile\n");
		exit(360);
	}
	SMspp		SP;
	stringPM	A(&SP);

	smsprun *run;
	run = NULL;

	StringmolSpatialConfigureFromFile(argv[2],&A,&run,1);

	s_ag *ag;

	for(ag=A.nowhead; ag != NULL; ag = ag->next){
		if(ag->status == B_ACTIVE){
			int found =0;
			comm_node *pb,*pend;
			pend=bbb;
			for(pb = bbb; pb != NULL; pb = pb->next){
				if(ag->spp->spp == pb->aspno){
					if(ag->pass->spp->spp == pb->pspno){
						found = 1;
						pb->count++;
						break;
					}
				}
				pend=pb;
			}
			if(!found){
				if(bbb == NULL)
					bbb = ReactionNewObservedPairing(ag->spp->spp,ag->pass->spp->spp,ag->spp->S,ag->pass->spp->S);
				else
				pend->next = ReactionNewObservedPairing(ag->spp->spp,ag->pass->spp->spp,ag->spp->S,ag->pass->spp->S);
			}
		}
	}

	comm_node *pb;
	FILE *fp;
	fp = fopen("1mcommunity.dat","w");

	printf("ACT\tPASS\tCOUNT\tSEQA\tSEQB\n");
	fprintf(fp,"ACT,PASS,COUNT,SEQA,SEQB\n");
	for(pb = bbb; pb != NULL; pb = pb->next){
		printf("%d\t%d\t%d\t%s\t%s\n",pb->aspno,pb->pspno,pb->count,pb->Sa,pb->Sp);
		fprintf(fp,"%d,%d,%d,%s,%s\n",pb->aspno,pb->pspno,pb->count,pb->Sa,pb->Sp);
	}
	printf("Finished!\n");
	fflush(stdout);
	fclose(fp);
	return 1;
}
