/* Copyright (C) 2009-2012 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#include "trigutil.h"
extern "C" {

	#include "memoryutil.h" //for memerror
	#include "randutil.h"

	#include "params.h"
}

#include "rules.h"
#include "agents_base.h"





/******************************************************************************
* @brief calculate the chances of two agents being close enough to bind
*
* @details propensity is aspatial-stringmol's way of calculating the likelihood
*          of two molecules meeting - this function implements the equation
*
* @param[in] n number of agents
*
* @return 1 if reaction can happen; 0 if not
*****************************************************************************/
int agents_base::PropensityEquation(const int n){

	float agarea = (float) M_PI*pow((float) agrad,2);
	float cellarea = M_PI*pow((float) vcellrad-(2.*(float) move),2);
	float arearatio = agarea/cellarea;


	if(n){
		//reactant coverage = 1 - ( 1- (area of reactant/area of cell) )^(number of reactants)
		float cov = 1.-pow(1.-arearatio,n) ;

		//printf("%d\t%f\t%f\t%f\t%f\n",n,cov,agarea,cellarea,arearatio);

		float rno = RandomBetween0And1();

		if(rno<cov)
			return 1;
		else
			return 0;
	}
	else{
		//printf("Zero mean coverage\n");
		return 0;
	}
}





//creators and destructors
agents_base::agents_base(){

	//Need to nullify these because stringPM doesn't use them
	//aat=NULL;
	//aac=NULL;
	//adc=NULL;
	//aro=NULL;

	//for preset
	bmax = 100;
	bct = (int *) malloc(bmax * sizeof(int));
	bpp = (int *) malloc(bmax * sizeof(int));
	memset(bct,0,bmax*sizeof(int));
	memset(bpp,0,bmax*sizeof(int));
	fr= NULL;
	tr = NULL;
	
	nsteps = 0;

	/*
	int i,j,k,found;
	float x1,y1,x2,y2;
	float dist,rad = 400-(9.9*2);
	for(i=0;i<bmax;i++){
		printf("Doing %d\n",i);
		for(j=0;j<100000;j++){
			if(!j%10000){
				printf("*");fflush(stdout);
			}
			found=0;
			rand_in_rad(rad,&x1,&y1);
			for(k=0;k<i;k++){
				rand_in_rad(rad,&x2,&y2);
				dist = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
				if(dist<9.9){
					found = 1;
				}
			}
			note_propensity(i,found);
		}
		printf("\n");
	}

	print_propensity(stdout);
	exit(7);
	*/

	preset();

}

agents_base::~agents_base(){

	clearout(0);

	free(bct);
	free(bpp);
}



/******************************************************************************
* @brief keep a record of this propensity event
*
* @details propensity is aspatial-stringmol's way of calculaing the likelihood
*          of two molecules meeting
*
* @param[in] N the number todo: what does it do?
*
* @param[in] X another number todo: what does it do?
*****************************************************************************/
void agents_base::PropensityRecord(int N,int X){
	if(N<bmax){
		bct[N]++;
		if(X)
			bpp[N]++;
	}

}


void agents_base::print_propensity(FILE *fp){
	int i;
	fprintf(fp,"\nPropensity table\n");
	for(i=0;i<bmax;i++)
		fprintf(fp,"%d\t%d\t%d\n",i,bct[i],bpp[i]);
	fflush(fp);
}


void agents_base::preset(){

	ifxhead = NULL;
	//btab = NULL;
	//com = NULL;
	//dcom = NULL;


	cellrad = 2500;
	agrad = 10;
	move = 0;
	energy = 0;
	//divtime = 0;
	vcellrad = 0;
}




// VJH - added this function for youShare - some alterations on merge
void agents_base::ConfigLoad(const char *fn, char *fninput, int test=0, int verbose=0){

	ParametersLoad(fn,test,verbose);

	// VJH load_agents(fn,test);
	AgentsLoad(fn,fninput,test);
}




/******************************************************************************
* @brief Load the parameters set from a config file
*
* @details loads the run parameters
*
* @param[in] fn name of the config file
*
* @param[in] test poorly designed test flag - todo: remove
*
* @param[in] verbose flag for verbose output
*
* @return 0 if error; 1 if success... todo reverse this
*****************************************************************************/
int agents_base::ParametersLoad(const char *fn, int test, int verbose){

	FILE *fp;

	if((fp=fopen(fn,"r"))!=NULL){
		float tmpen;
		int err = 0;
		int e=  ParameterReadFloat(fp,"CELLRAD",&cellrad, verbose);
		if(e>1)err++;

		vcellrad = cellrad;
		e=  ParameterReadFloat(fp,"AGRAD",&agrad, verbose);
		if(e>1)err++;

		//err +=  ParameterReadFloat(fp,"MOVE",&move);
		e=  ParameterReadFloat(fp,"ENERGY",&tmpen, verbose);
		if(e>1)err++;
		energy = (int) tmpen;

		e=  ParameterReadFloat(fp,"NSTEPS",&nsteps, verbose);
		if(e>1)err++;
		//err += 	ParameterReadFloat(fp,"DIVTIME",&divtime);

		if(err){
			printf("Some error reading config file\n");
			fclose(fp);
			return 0;
		}
		else{
			if(test){//Make sure things will collide
				vcellrad=cellrad=agrad;
			}
			fclose(fp);
			return 1;
		}
	}
	else{
		printf("Unable to open file %s\n",fn);
		fflush(stdout);
		return 0;
	}

}







void agents_base::clearout(int verbose){

	s_ix	*ixp,*ixp2;

	if(verbose){
		printf("Starting agents_base clearout..");fflush(stdout);
	}
	
	ixp=ifxhead;
	while(ixp!=NULL){
		ixp2=ixp->next;
		free(ixp);
		ixp=ixp2;
	}

	preset();

	if(verbose){printf("....finished\n");fflush(stdout);}
}
