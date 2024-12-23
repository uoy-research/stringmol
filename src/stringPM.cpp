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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "memoryutil.h"
#include "mt19937-2.h"
#include "randutil.h"
#include "params.h"
#include "hsort.h"

//string stuff
#include "stringmanip.h"
#include "alignment.h"
#include "instructions.h"

//metabolism stuff
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"

#include "error_codes.h"

//For debugging:
//#define VERBOSE
//#define V_VERBOSE

//Sticky stringmol:
//#define BIND_ALL

//For decay experiments in stringmol ALife paper:
//#define LONG_DECAY
//#define UNB_DECAY_ONLY

//Brutal hack for line length
//extern const int maxl = 2000;
//extern const int maxl0 = maxl+1; //allow room for a terminating 0




stringPM::stringPM(SMspp * pSP){

	if(pSP!=NULL)
		spl = pSP;
	else{
		spl = new SMspp();
	}

	dodecay=1;
	loadtype = L_REPLICABLE;

	grid = NULL;

	swlist=NULL;

	blosum = NULL;
	blosum = (swt *) malloc(sizeof(swt));
	blosum->N=0;
	blosum->T=NULL;
	blosum->key=NULL;
	SetHeadsAndDefaults();
	agct = 0;
	//species=NULL;
	spp_count=1;
	verbose_bind=0;

	//Set defaults:

	maxl = 2000;
	maxl0 = maxl+1; //allow room for a terminating 0
	estep = 20;

	/** This is a toggle to turn the '+' operator on and off
	 *
	 * It can be set using the 'GRANULAR' flag in the config file
	 *
	 */
	granular_1=0;

	splprint = 10000;

	report_every = 10000;	//How often to write splists and configs
	image_every	 = 100; //How often to generate an image (spatial stringmol only)


}


stringPM& stringPM::operator=(const stringPM &spm){

	spl      = spm.spl;

	dodecay  = spm.dodecay;
	loadtype = spm.loadtype;

	grid     = spm.grid;

	swlist	 = spm.swlist;

	blosum = NULL;
	blosum = (swt *) malloc(sizeof(swt));
	blosum->N   = 0;
	blosum->T   = NULL;
	blosum->key = NULL;
	
	SetHeadsAndDefaults();
	
	agct = spm.agct;
	spp_count = spm.agct;
	verbose_bind = spm.verbose_bind;

	//Set defaults:

	maxl = spm.maxl;
	maxl0 = spm.maxl0; //allow room for a terminating 0
	estep = spm.estep;

	/** This is a toggle to turn the '+' operator on and off
	 *
	 * It can be set using the 'GRANULAR' flag in the config file
	 *
	 */
	granular_1=spm.granular_1;

	splprint = spm.splprint;

	report_every = spm.report_every;   //How often to write splists and configs
	image_every  = spm.image_every;	   //How often to generate an image (spatial stringmol only)


	return *this;
}


//Copy constructor uses operator=
stringPM::stringPM(const stringPM& spm){
	*this = spm;
}





stringPM::~stringPM() {
	BucketReset(0);
}




//TODO(sjh): refactor this properly!
/*******************************************************************************
* @brief prepare structs and defaults
*******************************************************************************/
void stringPM::SetHeadsAndDefaults(){

	nowhead = NULL;
	nexthead = NULL;

	//note - it's possible that this'll be called twice - but will do no harm!
	agents_base::ParametersSetDefaults();
}





/* Make sure the masseage memory is freed - so you can't do printf("%s",parse_error(1)); - the memory will not be freed.
 *
 */
char * stringPM::parse_error(int errno){

	char * message;
	message = (char *) malloc(100*sizeof(char));

	memset(message,0,100*sizeof(char));

	sprintf(message,"Unspecified error");

	return message;

}





float stringPM::load_mut(const char *fn, int verbose){

	FILE *fp;
	float mut;

	if((fp=fopen(fn,"r"))!=NULL){

		char *emsg;
		int finderr=ParameterReadFloat(fp,"MUTATE",&mut,verbose);
		fclose(fp);

		switch(finderr){
		case 0:
			subrate = mut;
			indelrate = mut;
			if(mut<FLT_MIN)
				domut=0;
			else
				domut=1;
			if(verbose)printf("MUTATE setting: subrate = %f, indelrate = %f\n",subrate,indelrate);
			return 0;
			break;
		case 1:
			printf("No mutation rate found in config file\nUsing ALifeXII values instead\n");
			indelrate = 0.0000000306125;
			subrate=0.00001;
			domut=1;
			return 0;
			break;
		default: //Some other error
			printf("Error %d(%s) in loading mutation scheme\n",finderr,emsg=parse_error(finderr));
			free(emsg);
			return 1;
			break;

		}
	}
	return 2;
}





float stringPM::load_decay(const char *fn, int verbose){

	FILE *fp;
	float dec;

	if((fp=fopen(fn,"r"))!=NULL){

		int finderr=ParameterReadFloat(fp,"DECAY",&dec,verbose);
		fclose(fp);

		switch(finderr){
		case 0:
			decayrate = dec;
			if(verbose)printf("DECAY rate setting = %f\n",decayrate);
			return 0;
		case 1:
			printf("No decay rate found in config file\nUsing ALifeXII values instead\n");
			decayrate = 1.0/pow(65,2);
			return 0;
			break;
		default: //Some other error
			printf("Error %d in loading decay scheme\n",finderr);
			return 1;
			break;

		}
	}
	return 2;
}




//TODO: This should be in alignment.cpp
int stringPM::load_table_matrix(const char *fn){
	FILE *fp;
	printf("File name is %s",fn);
	if((fp=fopen(fn,"r"))!=NULL){
		int i,j; 
		char *p;
		char line[maxl];
		char label[maxl];
		
		if(fgets(line,maxl,fp)==NULL){
			printf("ERROR in reading line during load_table_matrix() in stringPM.cpp");
		}
		memset(label,0,maxl);
		sscanf(line,"%2000s",label);

		//Set N and key:
		blosum->N = strlen(label);
		blosum->key = (char *)malloc((blosum->N+1) * sizeof(char));
		memset(blosum->key,0,(blosum->N+1)*sizeof(char));
		strncpy(blosum->key,label,blosum->N);

		//Set up the data structure
		blosum->T = (float **) malloc((blosum->N +1) * sizeof(float *));
		for(i=0;i<blosum->N+1;i++){
			blosum->T[i] = (float *) malloc(blosum->N * sizeof(float));
		}

		for(i=0;i<blosum->N+1;i++){
			if(fgets(line,maxl,fp)==NULL){
				printf("ERROR in reading line during load_table_matrix() in stringPM.cpp");
			}
			p=strtok(line,", \t");
			for(j=0;j<blosum->N;j++){
				sscanf(p,"%f",&(blosum->T[i][j]));
				p = strtok(NULL,", \t");
			}
		}

		blosum->adj = (int **) malloc((blosum->N) * sizeof(int *));
		for(i=0;i<blosum->N;i++){
			blosum->adj[i] = (int *) malloc(blosum->N * sizeof(int));
		}
		for(i=0;i<blosum->N;i++){
			if(fgets(line,maxl,fp)==NULL){
				printf("ERROR in reading line during load_table_matrix() in stringPM.cpp\n");
			}
			p=strtok(line,", \t");
			for(j=0;j<blosum->N;j++){
				sscanf(p,"%d",&(blosum->adj[i][j]));
				p = strtok(NULL,", \t");
			}
		}

		fclose(fp);
		return 0;
	}
	else
		return 60;
}


int stringPM::load_table(const char *fn){


	//const int maxl=256;
	FILE *fp;


	if((fp=fopen(fn,"r"))!=NULL){
		int found = 0;
		char line[maxl];
		char label[maxl];
		char fn2[maxl];
		FILE *fp2;
		
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%2000s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"USING",5)){
				sscanf(line,"%*s %2000s",fn2);
				strcpy(swt_fn,fn2);
				found = 1;
				break;
			}
			if(!strncmp(line,"SUBMAT",6)){
				sscanf(line,"%*s %2000s",fn2);
				strcpy(swt_fn,fn2);
				found = 2;
				break;
			}
		}
		fclose(fp);
		if(!found){
			printf("No MIS or MTX file name found - using default values\n");
			sprintf(swt_fn,"UsingDefaults!");
			found = 3;
		}


		switch(found){
		case 1:
			if((fp2=fopen(fn2,"r"))!=NULL){
				int found2 = 0;
				//BUG: For some reason, when we open this again, we aren't at the beginning of the file..
				//can't understand why this is suddenly happening, but let's put a rewind in to fix it for now!
				rewind(fp2);
				while((fgets(line,maxl,fp2))!=NULL){
					memset(label,0,maxl);
					sscanf(line,"%2000s",label);
					//printf("line = %s",line);
					if(!strncmp(line,"SET",3)){
						sscanf(line,"%*s %2000s",fn2); //using fn2 to temporarily hold the alphabet...
						blosum->N = strlen(fn2);
						blosum->key = (char *)malloc((blosum->N+1) * sizeof(char));
						memset(blosum->key,0,(blosum->N+1)*sizeof(char));
						strncpy(blosum->key,fn2,blosum->N);

						blosum->T = (float **) malloc((blosum->N +1) * sizeof(float *));
						for(int i=0;i<blosum->N+1;i++){
							blosum->T[i] = (float *) malloc(blosum->N * sizeof(float));
						}

						table_from_string(blosum->T,blosum->key,blosum->N);

						//found2 = 1;
						fclose(fp2);

						//SUCCESS
						return 0;
					}
				}
				if(!found2){
					fclose(fp2);
					printf("No MIS string found in MIS file\n");
					return 2;
				}

			}
			else{
				printf("Unable to open MIS file %s",fn2);
				return 3;
			}
			break;
		case 2:
			return load_table_matrix(swt_fn);
			break;
		case 3:
			blosum =  default_table();
			return 0;
			break;
		}
		print_swt(stdout,blosum);
	}
	else
		return 4;

	return 5;


}


int stringPM::load_splist(const char *fn,int verbose){

	const int llen = maxl0 + 256;
	FILE *fp;
	char line[llen];
	char label[llen];
	int nspp=0;
	if((fp=fopen(fn,"r"))!=NULL){
		while((fgets(line,llen,fp))!=NULL){
			memset(label,0,llen);
			sscanf(line,"%2000s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"SPECIES",6)){
				nspp++;
			}
		}
	}

	int spno;

	if(nspp){//We've got species to create
		rewind(fp);

		while((fgets(line,llen,fp))!=NULL){
			memset(label,0,llen);
			sscanf(line,"%2000s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"SPECIES",6)){

				//Scan the line for sequence and number...
				memset(label,0,llen);
				sscanf(line,"%*s %d %2000s",&spno,label);
				//TODO: make sure extit is set at this point
				spl->getspp_from_string(label,timestep,maxl0,spno);
			}
		}
		fclose(fp);
	}
	else{
		fclose(fp);
	}
	return 0;
}


int stringPM::id_spp(l_spp *sp, s_ag *pag, int  aspno, char *spp_string){
	if(sp->spp == aspno){
		if(!strncmp(sp->S,spp_string,strlen(spp_string))){
			pag->spp = sp;
		}
		else{
			printf("ERROR species list is corrupted - numbers don't match strings\n");
			exit(0);
		}
	}
	return 0;
}



int stringPM::load_reactions(const char *fn, char *fntab, int test, int verbose){


	const int llen = maxl0 + 256;
	FILE *fp;
	char line[llen];

	s_ag *pag,*bag;

	int aspno,pspno,ano,pno;
	int it,iap,ipp,rt,rap,rpp,wt,wap,wpp,ft,fap,fpp, gx, gy;

	//First, let's see what format the reactions are in - look for the keyword "NUMAGENTS"


	if((fp=fopen(fn,"r"))!=NULL){
		while((fgets(line,llen,fp))!=NULL){

		}
		fclose(fp);
	}



	if((fp=fopen(fn,"r"))!=NULL){
		int linecount = 0;
		while((fgets(line,llen,fp))!=NULL){
			linecount++;
			
			char active_spp_string[llen];
			char passive_spp_string[llen];
			char active_string[llen];
			char passive_string[llen];
	
			memset(active_spp_string,0,llen);
			memset(active_string,0,llen);
			memset(passive_string,0,llen);
			memset(passive_spp_string,0,llen);
			sscanf(line,"%2000s",active_spp_string);
			//printf("line = %s",line);
			if(!strncmp(line,"REACTION",8)){
				/* We have a multi-line definition to deal with now... hold tight... */

				/* First line is the species of the active molecule */
				if((fgets(line,llen,fp))!=NULL){
					linecount++;
					sscanf(line,"%*s %d %2000s",&aspno, active_spp_string);
					//TODO: check that the species number and sequence are identical
				}
				else{
					printf("ERROR READING REACTION (active species) at line %d\n",linecount+1);
					return 1;
				}

				/* Second line is the state of the active molecule */
				if((fgets(line,llen,fp))!=NULL){
					linecount++;
					gx = gy = -1;
					sscanf(line,"%*s %d %2000s irwf: %d %d %d %d %d %d %d %d %d %d %d %d grid: %d %d",&ano, active_string,
							&it,&iap,&ipp,&rt,&rap,&rpp,&wt,&wap,&wpp,&ft,&fap,&fpp,
							&gx,&gy);
				}
				else{
					printf("ERROR READING REACTION (active state) at line %d\n",linecount+1);
					return 1;
				}

				pag = AgentMake('X');//TODO: fix this need for ascii codes... we have species numbers now!
				pag->S =(char *) malloc(maxl0*sizeof(char));
				memset(pag->S,0,maxl0*sizeof(char));
				strncpy(pag->S,active_string,maxl);//active_string));
				pag->len = strlen(pag->S);

				if(grid){
					pag->x = gx;
					pag->y = gy;
				}

				/* Third line is the state of the passive molecule */
				if((fgets(line,llen,fp))!=NULL){
					linecount++;
					gx = gy = -1;
					sscanf(line,"%*s %d %2000s grid: %d %d",&pno, passive_string, &gx, &gy);
				}
				else{
					printf("ERROR READING REACTION (passive state) at line %d\n",linecount+1);
					return 1;
				}
				bag = AgentMake('X');//TODO: fix this need for ascii codes... we have species numbers now!
				bag->S =(char *) malloc(maxl0*sizeof(char));
				memset(bag->S,0,maxl0*sizeof(char));
				strncpy(bag->S,passive_string,maxl);//passive_string));
				bag->len = strlen(bag->S);

				if(grid){
					bag->x = gx;
					bag->y = gy;
				}

				/* Fourth line is the species of the passive molecule */
				if((fgets(line,llen,fp))!=NULL){
					linecount++;
					sscanf(line,"%*s %d %2000s",&pspno, passive_spp_string);
					//TODO: check that the species number and sequence are identical
				}
				else{
					printf("ERROR READING REACTION (active species) at line %d\n",linecount+1);
					return 1;
				}


				// set the states - there's a function for this called set_exec, but we can't use that because the states may now

				pag->status = B_ACTIVE;
				pag->exec = NULL;
				pag->pass = bag;

				bag->status = B_PASSIVE;
				bag->exec = pag;
				bag->pass = NULL;

				//  Now we can set the pointers:
				pag->it = it;
				pag->i[0] = &(bag->S[iap]);
				pag->i[1] = &(pag->S[ipp]);
				pag->rt = rt;
				pag->r[0] = &(bag->S[rap]);
				pag->r[1] = &(pag->S[rpp]);
				pag->wt = wt;
				pag->w[0] = &(bag->S[wap]);
				pag->w[1] = &(pag->S[wpp]);
				pag->ft = ft;
				pag->f[0] = &(bag->S[fap]);
				pag->f[1] = &(pag->S[fpp]);

				//set the species data:
				l_spp *sp;
				for(sp = spl->species_list; sp!=NULL; sp=sp->next){
					id_spp(sp, pag, aspno,  active_spp_string);
					id_spp(sp, bag, pspno, passive_spp_string);
				}

				//Finally, append the reacting molecules to the simulation
				AgentAppend(&nowhead,pag);
				AgentAppend(&nowhead,bag);
			}
		}
	}
	return 0;
}






//TODO: check that the file pointer can be passed in like this!
s_ag * stringPM::read_unbound_agent(FILE **fp, char line[], const int llen){

	s_ag * pag;
	char code;
	int nag;
	char label[llen];


	memset(label,0,llen);
	sscanf(line,"%*s %2000s %d %c",label,&nag,&code);

	//bool getgridinfo = true;

	//make the agent
	//for(i=0;i<nag;i++){
	l_spp *s;


	pag = AgentMake(code);//,1);

	pag->S =(char *) malloc(maxl0*sizeof(char));

	memset(pag->S,0,maxl0*sizeof(char));

	if(strlen(label) > maxl0){
		printf("Unable to allocate enough space for this agent: \n%s\nConsider specifying MAXL in your config\n",label);
	}

	strncpy(pag->S,label,strlen(label));
	pag->len = strlen(pag->S);

	/*Load the coordinates on the grid for this agent
	 * THESE ENTRIES MUST BE IMMEDIATELY AFTER THE AGENT
	 *
	 *
	 * */
	if(grid){// && getgridinfo){
		if((fgets(line,llen,*fp))!=NULL){//sscanf(line,"%s",label);
			if(!strncmp(line,"GRIDPOS",7)){
				sscanf(line,"%*s %d %d ",&(pag->x),&(pag->y));
			}
			//else{
			//	getgridinfo = false;
			//}
		}
	}


	//No parents for these initial agents!
	pag->pp = spl->ParentsMake(NULL,NULL);

	//if(!i){

	//int stringPM::SpeciesListUpdate(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp)
	SpeciesListUpdate(pag,'I',1,NULL,NULL,0);
	s = spl->getspp(pag,timestep,maxl0);
	//TODO: tidy up handling of seed species, but for now:
	s->tspp = 0;
	//}

	//Record that the replicase copies itself
	pag->spp=s; //TODO: is this true for non-replicase agents..?

	return pag;

}






s_ag * stringPM::read_active_agent(FILE **fp, char line[], const int llen, int &pidx){

	int aspno;
	char active_spp_string[llen];
	char active_string[llen];
	int ano,
		it,iap,ipp,
		rt,rap,rpp,
		wt,wap,wpp,
		ft,fap,fpp,
		gx,gy;

	s_ag *pag;

	/* First line is the species of the active molecule */
	if((fgets(line,llen,*fp))!=NULL){
		linecount++;
		sscanf(line,"%*s %d %2000s",&aspno, active_spp_string);
		//TODO: check that the species number and sequence are identical
	}
	else{
		printf("ERROR READING REACTION (active species) at line %d\n",linecount+1);
		return NULL;
	}

	/* Second line is the state of the active molecule */
	if((fgets(line,llen,*fp))!=NULL){
		linecount++;
		gx = gy = -1;
		sscanf(line,"%*s %d %2000s irwf: %d %d %d %d %d %d %d %d %d %d %d %d grid: %d %d  pass_index: %d",
				&ano, active_string,
				&it,&iap,&ipp,&rt,&rap,&rpp,&wt,&wap,&wpp,&ft,&fap,&fpp,
				&gx,&gy,
				&pidx);
	}
	else{
		printf("ERROR READING REACTION (active state) at line %d\n",linecount+1);
		return NULL;
	}

	pag = AgentMake('X');//TODO: fix this need for ascii codes... we have species numbers now!
	pag->S =(char *) malloc(maxl0*sizeof(char));
	memset(pag->S,0,maxl0*sizeof(char));
	strncpy(pag->S,active_string,maxl);//active_string));
	pag->len = strlen(pag->S);

	if(grid){
		pag->x = gx;
		pag->y = gy;
	}

	/* Final step at this point is to set the pointers, but
	 * since we don't have the passive agent yet, we'll set the passive
	 * pointers to the equivalent position on the active string, and
	 * move it during the final stage of the load.
	 */

	//  Now we can set the pointers:
	pag->it = it;
	pag->i[0] = &(pag->S[iap]);
	pag->i[1] = &(pag->S[ipp]);
	pag->rt = rt;
	pag->r[0] = &(pag->S[rap]);
	pag->r[1] = &(pag->S[rpp]);
	pag->wt = wt;
	pag->w[0] = &(pag->S[wap]);
	pag->w[1] = &(pag->S[wpp]);
	pag->ft = ft;
	pag->f[0] = &(pag->S[fap]);
	pag->f[1] = &(pag->S[fpp]);

	//set the species data:
	l_spp *sp;
	for(sp = spl->species_list; sp!=NULL; sp=sp->next){
		id_spp(sp, pag, aspno,  active_spp_string);
		//id_spp(sp, bag, pspno, passive_spp_string);
	}

	return pag;
}





s_ag * stringPM::read_passive_agent(FILE **fp, char line[], const int llen){

	int gx,gy;
	int pno,pspno;
	char passive_spp_string[llen];
	char passive_string[llen];

	/* Third line is the state of the passive molecule */
	if((fgets(line,llen,*fp))!=NULL){
		linecount++;
		gx = gy = -1;
		sscanf(line,"%*s %d %2000s grid: %d %d",&pno, passive_string, &gx, &gy);
	}
	else{
		printf("ERROR READING REACTION (passive state) at line %d\n",linecount+1);
		return NULL;
	}

	/* Fourth line is the species of the passive molecule */
	if((fgets(line,llen,*fp))!=NULL){
		linecount++;
		sscanf(line,"%*s %d %2000s",&pspno, passive_spp_string);
		//TODO: check that the species number and sequence are identical
	}
	else{
		printf("ERROR READING REACTION (passive species) at line %d\n",linecount+1);
		return NULL;
	}

	s_ag *bag;

	/* Now create the molecule */
	bag = AgentMake('X');//TODO: fix this need for ascii codes... we have species numbers now!
	bag->S =(char *) malloc(maxl0*sizeof(char));
	memset(bag->S,0,maxl0*sizeof(char));
	strncpy(bag->S,passive_string,maxl);//passive_string));
	bag->len = strlen(bag->S);

	if(grid){
		bag->x = gx;
		bag->y = gy;
	}

	//set the species data:
	l_spp *sp;
	for(sp = spl->species_list; sp!=NULL; sp=sp->next){
		//id_spp(sp, pag, aspno,  active_spp_string);
		id_spp(sp, bag, pspno, passive_spp_string);
	}

	return bag;
}





smsprun * stringPM::init_smprun(const int gridx, const int gridy){

	grid = (smsprun *) malloc(sizeof(smsprun));

	//TODO: make grid size changeable via config...
	grid->gridx=gridx;
	grid->gridy=gridy;

	grid->grid=(s_ag ***) malloc(grid->gridx*sizeof(s_ag **));
	grid->status=(s_gstatus **) malloc(grid->gridx*sizeof(s_gstatus *));

	for(int i=0;i<grid->gridx;i++){
		grid->grid[i] = (s_ag **) malloc(grid->gridy*sizeof(s_ag *));
		//grid->status[i] = (s_gstatus *) malloc(grid->gridy*sizeof(s_gstatus));
		grid->status[i] = static_cast<s_gstatus *> (malloc(grid->gridy*sizeof(s_gstatus)));

		for(int j=0;j<grid->gridy;j++){
			grid->grid[i][j]=NULL;
			grid->status[i][j]=G_EMPTY;
		}
	}

	return grid;
}





int stringPM::load_replicable(const char *fn){
	int nag=0;
	const int llen = maxl0 + 256;
	FILE *fp;
	s_ag **agarray;
	int *parray;
	
	agarray = NULL;
	parray = NULL;

	if((fp = fopen(fn,"r"))!=NULL){

		char line[llen];
		char label[llen];
		while((fgets(line,llen,fp))!=NULL){
			memset(label,0,llen);
			sscanf(line,"%2000s",label);
			if(!strncmp(line,"NUMAGENTS",9)){
				sscanf(line,"NUMAGENTS %d",&nag);
				break;
			}
		}

		if(nag < 1){
			printf("Number of agents not specified - impossible to configure using repclicable method...\n");
			printf("Reverting to original load method (pre 2017)...\n");
			return(234);
		}

		//Set up a 'scaffold' array from which we'll build the linked list...
		agarray = (s_ag **) malloc(nag * sizeof( s_ag *));
		parray = (int *) malloc(nag * sizeof(int));
		memset(parray,0,nag*sizeof(int));

		for(int ii =0; ii< nag; ii++){
			agarray[ii] = NULL;
		}

		//we *have* to assume that the agents are in order.
		int seenag = 0;

		while((fgets(line,llen,fp))!=NULL){
			memset(label,0,llen);
			sscanf(line,"%2000s",label);

			//Work out the bind state - different load methods for each.
			s_bind bs;
			s_ag *pag;

			bool found = false;
			if(!strncmp(label,"AGENT",5)){
				bs = B_UNBOUND;
				found = true;
			}
			if(!strncmp(label,"REACTION",8)){
				/*TODO: this is a bit risky because there are lines beginning REACTION_ACTIVE
				 * and REACTION_PASSIVE.. should be ok if we process correctly though... */
				bs = B_ACTIVE;
				found = true;
			}
			if(!strncmp(label,"PASSIVE",7)){
				bs = B_PASSIVE;
				found = true;
			}

			if(found){
				int pidx=0;
				switch(bs){
				case B_UNBOUND:
					pag = read_unbound_agent(&fp, line,llen);
					pag->status = B_UNBOUND;
					break;
				case B_ACTIVE:
					pag = read_active_agent(&fp, line, llen, pidx);
					pag->status = B_ACTIVE;
					parray[seenag] = pidx;
					break;
				case B_PASSIVE:
					pag = read_passive_agent(&fp, line,llen);
					pag->status = B_PASSIVE;
					break;
				}
				agarray[seenag] = pag;
				seenag++;
			}

		}
	fclose(fp);
	}

	s_ag *bag, *pag;
	int pp;

	//TODO: Now we need to go through and get all the active mols to point at the passive ones.
	for(int ii=0;ii<nag;ii++){

		//add the agent to the list...
		AgentAppend(&nowhead,agarray[ii]);

		if(agarray[ii]->status == B_ACTIVE){

			pag = agarray[ii];

			//Get the index of the passive partner..
			pag->pass = agarray[parray[ii]];
			bag = pag->pass;
			pag->pass->exec = pag;

			//  Now we can set the pointers (adapted from load_reactions)
			//agarray[ii]->it = it;
			pp = pag->i[0] - pag->S;
			pag->i[0] = &(bag->S[pp]);

			//agarray[ii]->i[1] = &(agarray[ii]->S[ipp]);
			//agarray[ii]->rt = rt;
			pp = pag->r[0] - pag->S;
			pag->r[0] = &(bag->S[pp]);
			//agarray[ii]->r[1] = &(agarray[ii]->S[rpp]);
			//agarray[ii]->wt = wt;


			pp = pag->w[0] - pag->S;
			pag->w[0] = &(bag->S[pp]);
			//agarray[ii]->w[1] = &(agarray[ii]->S[wpp]);
			//agarray[ii]->ft = ft;


			pp = pag->f[0] - pag->S;
			pag->f[0] = &(bag->S[pp]);
			//agarray[ii]->f[1] = &(agarray[ii]->S[fpp]);

		}
	}

	free(agarray);
	free(parray);

	return 1;

}





/******************************************************************************
* @brief Load the molecule set from a config file
*
* @details loads the agents, but *also* many of the run parameters (todo: fix this).
*          exit on failure to handle essential parameters (todo: proper error codes!)
*
* @param[in] fn name of the config file
*
* @param[in] fntab name of the file containing the BLOSUM matrix
*
* @param[in] test poorly designed test flag - todo: remove
*
* @param[in] verbose flag for verbose output
*
*
* @return 0 if error; 1 if success... todo reverse
*****************************************************************************/
int stringPM::AgentsLoad(const char *fn, char *fntab, int test, int verbose){


	//Load the maximum permitted length of the molecules (or the default)
	int dmxl = maxl;
	int err = ParameterReadOrDefineUnsignedInt(fn,"MAXLEN", &maxl, dmxl, 0);
	if(err>1){
		printf("ERROR %d on loading max line length (MAXLEN)\n",err);
		exit(0);
	}
	maxl0 = maxl+1;

	timestep = 0;
	err = ParameterReadOrDefineUnsignedInt(fn,"EXTIT", &timestep, 0, 0);
	if(err>1){
		printf("ERROR %d on loading external iteration value (EXTIT)\n",err);
		exit(PARAM_LOAD_ERROR);
	}

	err = ParameterReadOrDefineUnsignedInt(fn,"REPORTEVERY", &report_every, report_every, 0);
	if(err>1){
		printf("ERROR %d on loading reporting frequency\n",err);
		exit(PARAM_LOAD_ERROR);
	}


	err = ParameterReadOrDefineUnsignedInt(fn,"IMAGEEVERY", &image_every, image_every, 0);
	if(err>1){
		printf("ERROR %d on loading imaging frequency\n",err);
		exit(PARAM_LOAD_ERROR);
	}



	//TODO: Not sure that this has ever been used - but it might be useful
	err = ParameterReadOrDefineUnsignedInt(fn,"SPLPRINT", &splprint, splprint, 0);
	if(err>1){
		printf("ERROR %d on loading specieslist print frequency\n",err);
		exit(PARAM_LOAD_ERROR);
	}
	maxl0 = maxl+1;


	const int llen = maxl0 + 256;

	char code;
	int nag;
	s_ag *pag;
	int ntt = 0;



	int destep = estep;
	int estep_err = ParameterReadOrDefineUnsignedInt(fn,"ESTEP", &estep, destep, 0);
	if(estep_err>1){
		printf("ERROR %d on loading energy per time step (ESTEP)\n",estep_err);
		exit(PARAM_LOAD_ERROR);
	}


	//ALSO important to load blosum table here too!
	int table_err;
	if(fntab)
		table_err = load_table_matrix(fntab);
	else
		table_err = load_table((const char *) fn);
	if(table_err){
		printf("ERROR %d on loading table\n",table_err);
	}

	int mut_err = load_mut(fn,verbose);
	if(mut_err){
		printf("ERROR %d on loading mutation rate\n",mut_err);
	}

	int decay_err = load_decay(fn,verbose);
	if(decay_err){
		printf("ERROR %d on loading decay rate\n",decay_err);
	}


	/*Load any grid info (if needed)
	 */

	int gdefault = 0;
	unsigned int gridx,gridy;
	err = ParameterReadOrDefineUnsignedInt(fn,"GRIDX", &gridx, gdefault, 0);
	if(err>1){
		printf("ERROR %d on loading gridx (GRIDX)\n",err);
		exit(PARAM_LOAD_ERROR);
	}
	else{
		if(gridx>0){//Then we've got a grid definition and we need a y
			err = ParameterReadOrDefineUnsignedInt(fn,"GRIDY", &gridy, gdefault, 0);
			if(err>1){
				printf("ERROR %d on loading gridy (GRIDY)\n",err);
				exit(PARAM_LOAD_ERROR);
			}
			if(gridy==0){
				printf("ERROR: GRIDX is %u, but GRIDY is 0\n",gridx);
				exit(PARAM_LOAD_ERROR);
			}
			init_smprun(gridx,gridy);
		}
	}


	ntt = 0;

	/* READ THE SPECIES LISTS (IF ANY)*
	 * Important to do this first, so that if there any agents that aren't on the list, it isn't cocked up
	 * Important also to maintaint the species numbering - will save a *lot* of hassle later!
	 */
	int spp_err = load_splist(fn,verbose);
	if(spp_err){//TODO: this will never be called because the function always returns 0
		printf("Problem loading species list...exiting\n");
		exit(456);
	}

	/* READ THE AGENTS IN NOW */

	/* See if we can use the replicable method; use the old method if not... */
	if(loadtype == L_REPLICABLE){
		if(load_replicable(fn) == 234)
			loadtype = L_EVOEVO;
	}

	/* This will only be called if we've reverted to the old method (or set it elsewhere) */
	if(loadtype == L_EVOEVO){/*Only used for evoevo 2016 runs: */

		//TODO: this is the actual agent-loading part of the function... rename / refactor needed...
		FILE *fp;
		
		if((fp=fopen(fn,"r"))!=NULL){
			char line[llen];
			char label[llen];
			while((fgets(line,llen,fp))!=NULL){
				memset(label,0,llen);
				//TODO: We need a better way of limiting the chars 
				//read into label - perhaps with a #define or two...
				sscanf(line,"%2000s",label);
				//printf("line = %s",line);
				if(!strncmp(line,"AGENT",5)){
					ntt++;
				}
			}

			rewind(fp);

			bool getgridinfo;

			while((fgets(line,llen,fp))!=NULL){

				memset(label,0,llen);
				sscanf(line,"%2000s",label);
				if(!strncmp(line,"AGENT",5)){

					memset(label,0,llen);
					sscanf(line,"%*s %2000s %d %c",label,&nag,&code);

					if(test){//todo -is this still necessary?
						nag=1;
					}

					getgridinfo = true;
						
						
					l_spp *s;
					s=NULL;
						

					//make the agents
					for(int i=0;i<nag;i++){


						pag = AgentMake(code);//,1);

						pag->S =(char *) malloc(maxl0*sizeof(char));

						memset(pag->S,0,maxl0*sizeof(char));

						if(strlen(label) > maxl0){
							printf("Unable to allocate enough space for this agent: \n%s\nConsider specifying MAXL in your config\n",label);
						}

						strncpy(pag->S,label,strlen(label));
						pag->len = strlen(pag->S);

						/*Load the coordinates on the grid for this agent
						 * THESE ENTRIES MUST BE IMMEDIATELY AFTER THE AGENT
						 *
						 *
						 * */
						if(grid && getgridinfo){
							if((fgets(line,llen,fp))!=NULL){//sscanf(line,"%s",label);
								if(!strncmp(line,"GRIDPOS",7)){
									sscanf(line,"%*s %d %d ",&(pag->x),&(pag->y));
								}
								else{
									getgridinfo = false;
								}
							}
						}


						//No parents for these initial agents!
						pag->pp = spl->ParentsMake(NULL,NULL);

						if(!i){

							//int stringPM::SpeciesListUpdate(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp)
							SpeciesListUpdate(pag,'I',1,NULL,NULL,0);
							s = spl->getspp(pag,timestep,maxl0);
							//TODO: tidy up handling of seed species, but for now:
							s->tspp = 0;
						}

						//Record that the replicase copies itself
						pag->spp=s; //TODO: is this true for non-replicase agents..?

						//printf("agent %d, %c, %0.3f %0.3f\n",i,pag->label,pag->x,pag->y);
						if(nowhead == NULL){
							nowhead = pag;
						}
						else{
							AgentAppend(&nowhead,pag);
						}
					}
				}
			}

			//close
			load_reactions(fn,fntab,test,verbose);
			fclose(fp);
			return 1;
		}
		else{
			printf("Unable to open file %s\n",fn);
			return 0;
		}

	}

	return 1;
}





/******************************************************************************
* @brief count species
*
* @return the number of species
*****************************************************************************/
int stringPM::SpeciesCount(){
	
	int nag, *done;
	s_ag *pa;
	int spc;
	nag=AgentsCount(nowhead,-1);
	int sppcount=0;

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	int finished;
	do{
		int i = 0;
		int found=0;
		
		finished = 1;
		for(i=0,pa=nowhead;i<nag;i++,pa=pa->next){
			if(!done[i]){
				if(!found){
					//Set the species
					spc = pa->spp->spp;
					//Check it off the list
					done[i]=1;
					sppcount++;
					finished=0;
					found=1;
				}
				else{
					if(pa->spp->spp==spc){
						done[i]=1;
					}
				}
			}
		}
	}while(!finished);
	free(done);
	return sppcount;
}



/******************************************************************************
* @brief print the number of species to file/screen
*
* @param[in] fp file pointer
*
* @param[in] style 0: print horizontal; 1: print vertical
*
* @param[in] state -1: print all, otherwise print only those agents whose
*            status == state
*****************************************************************************/
void stringPM::SpeciesPrintCount(FILE *fp,int style, int state){

	//char fn[128];
	s_ag *pa;
	int finished = 0;
	int nag,*done;
	int *spno, *spct, nspp = SpeciesCount();

	nag = AgentsCount(nowhead,-1);

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	spno = (int *) malloc(nspp*sizeof(int));
	memset(spno,0,nspp*sizeof(int));
	spct = (int *) malloc(nspp*sizeof(int));
	memset(spct,0,nspp*sizeof(int));
	int spidx=0;

	do{
		int i = 0;
		int found=0;
		finished = 1;
		for(i=0,pa=nowhead;i<nag;i++,pa=pa->next){
			if(state==-1 || pa->status == state){
				if(!done[i]){
					if(!found){
						done[i]=1;
						spno[spidx]=pa->spp->spp;
						spct[spidx]=1;
						finished=0;
						found=1;
					}
					else{
						if(pa->spp->spp==spno[spidx]){
							done[i]=1;
							spct[spidx]++;
						}
					}
				}
			}
		}
		spidx++;
	}while(!finished);

	int *spindx;
	spindx = (int *) malloc(nspp*sizeof(int));
	for(int ii=0;ii<nspp;ii++){
		spindx[ii]=ii;
	}

	if(nspp>1){

		//print_hsort_data(spno,spindx,nspp,stdout);
		HeapsortByIndex(spno,spindx,nspp);
	}
	else{
		spindx[0]=0;
	}

	switch(style){
	case 0: //Horizontal
		fprintf(fp,"-------");
		for(int i=0;i<nspp;i++){
			fprintf(fp,"-------");
		}
		fprintf(fp,"\n%d\t",(int) timestep);
		for(int i=0;i<nspp;i++){
			fprintf(fp,"%03d\t",spno[spindx[i]]);
		}
		fprintf(fp,"\n%d\t",(int) timestep);
		for(int i=0;i<nspp;i++){
			fprintf(fp,"%03d\t",spct[spindx[i]]);
		}
		fprintf(fp,"\n");
		break;
	case 1: //Vertical
		for(int i=0;i<nspp;i++){
			fprintf(fp,"%d,%d,%d\n",(int) timestep,spno[spindx[i]],spct[spindx[i]]);
		}
		break;
	}


	free(done);
	free(spno);
	free(spct);
	free(spindx);
}








/*TODO: Delete when we have found this functionality elsewhere in the code
void stringPM::get_spp_count(int state){

	int nag,*done;
	int *spno, *spct, nspp = count_spp();

	nag = nagents(nowhead,-1);

	l_spp *pSpp;

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	spno = (int *) malloc(nspp*sizeof(int));
	memset(spno,0,nspp*sizeof(int));
	spct = (int *) malloc(nspp*sizeof(int));
	memset(spct,0,nspp*sizeof(int));

	//Set the spl->count param to zero
	for(pSpp=spl->species;pSpp!=NULL;pSpp=pSpp->next){
		pSpp->count = 0;
	}


	s_ag *p;
	for(p = nowhead;p!=NULL;p=p->next){
		if(state==-1 || p->status == state){
			p->spp->count++;
		}
	}
}*/




/******************************************************************************
 * @brief append an agent to a list
 *
 * @param[in] list the list of agents
 *
 * @param[in] ag the agent
 *
 * @return 0 always
 *****************************************************************************/
int stringPM::AgentAppend(s_ag **list, s_ag *ag){
	s_ag *pag;

	//printf("appending: list = %p, ag = %p\n",*list,ag);
	if(*list==NULL){
		*list=ag;
		//printf("appended: list = %p, ag = %p\n",*list,ag);
	}
	else{
		pag = *list;
		while(pag->next != NULL){
			pag = pag->next;
		}
		pag->next = ag;
		ag->prev = pag;
	}

	return 0;
}





/*******************************************************************************
* @brief determine whether an agent is present in a list (by address)
*
* @param[in] list the list
*
* @param[in] tag the agent
*
* @return true or false
*******************************************************************************/
bool stringPM::AgentAddressInList(s_ag *list,const s_ag *tag){

	s_ag *pag;

	for(pag=list;pag!=NULL;pag=pag->next){
		if(tag == pag){
			return true;
		}
	}
	return false;
}




/******************************************************************************
 * @brief extract an agent from a list
 *
 * @param[in] list the list of agents
 *
 * @param[in] ag the agent
 *
 * @return 0 always
 *****************************************************************************/
int stringPM::AgentExtract(s_ag **list, s_ag *ag){

	//printf("extracting: list = %p, ag = %p\n",*list,ag);
	if(ag == *list){
		*list = ag->next;
		if(*list !=NULL)
			(*list)->prev = NULL;
	}
	else{
		if(ag->prev==NULL){
			printf("Error in extract_ag: No previous member of the list!\n");
		}else{	
			ag->prev->next = ag->next;
			if(ag->next != NULL)
				ag->next->prev = ag->prev;
		}
	}
	ag->prev=NULL;
	ag->next=NULL;
	return 0;
}





/******************************************************************************
 * @brief free memory used by an agent
 *
 * @param[in] pag the agent to free
 *
 * @return 0 always
 *****************************************************************************/
int stringPM::AgentFree(s_ag *pag){

	if(pag->S != NULL){
		//printf("destroying agent %d, code = %s\n",pag->idx,pag->S);
		free(pag->S);
	}

	free(pag);
	//TODO: we should set this to null at the moment we have to do it after each call to this function...!
	//pag = NULL;

	return 0;
}





/******************************************************************************
 * @brief count agents
 *
 * @param[in] head: usually 'nowhead' or 'nexthead'
 *
 * @param[in] state if -1: count all; else count with a particular status
 *
 * @return the count
 *****************************************************************************/
int stringPM::AgentsCount(s_ag *head, int state){
	s_ag *pag;
	int count=0;
	pag = head;
	while(pag!=NULL){
		switch(state){
		case -1:
			count++;
			break;
		default:
			if(pag->status == state)
				count++;
				/* no break */
		}
		pag=pag->next;
	}
	return count;
}





/*******************************************************************************
 * @brief select a random agent from a list
 *
 * @param[in] head: usually 'nowhead' or 'nexthead'
 *
 * @param[in] state if -1: count all; else count with a particular status
 *
 * @return an agent, or NULL if not available
 ******************************************************************************/
s_ag * stringPM::AgentSelectRandomly(s_ag *head, int state){
	int count = AgentsCount(head,state);
	int i,pos;
	s_ag *pag,**arr;

	if(!count)
		return NULL;

	pag=NULL;

	switch(state){
	case -1:
		while(pag==NULL){
			pos = (int) (count * RandomBetween0And1());
			//printf("count = %d, pos = %d\n",count,pos);
			pag = head;
			for(i=0;i<pos;i++){
				pag = pag->next;
			}
		}
		break;
	case B_UNBOUND:
	case B_ACTIVE:
	case B_PASSIVE:
		arr = (s_ag **) malloc(count*sizeof(s_ag *));
		i=0;
		pag=head;
		while(pag!=NULL){
			if(pag->status==state)
				arr[i++]=pag;
			pag=pag->next;
		}
		pos=count;
		while(pos==count){
			pos = (int) (count * RandomBetween0And1());
		}
		pag=arr[pos];
		free(arr);
		break;
	default:
		pag = NULL;

	}
	return pag;
}



/**
 * Create an 'agent', which is a string 'molecule'
 * NB: The string is not allocated here - done outside the function
 *
 * parameters:
 *
 * alab: a single-character label for the string (e.g. 'C')
 *
 */
/*******************************************************************************
 * @brief allocate memory and set up structure of a new agent
 *
 * @param[in] alab an integer identifier
 *
 * @return the agent
* *****************************************************************************/
s_ag * stringPM::AgentMake(int alab){

	s_ag *ag;

	//printf("Spatial make_ag called\n");fflush(stdout);

	//ATTENTION! this is how you convert from C-style to C++-style casts
	//https://embeddedartistry.com/blog/2017/03/15/c-casting-or-oh-no-they-broke-malloc/
	//if((ag = (s_ag *) mymalloc(1,sizeof(s_ag)))!=NULL){
	if(( ag = static_cast<s_ag *> (MallocOrExit(1, sizeof(s_ag))))!=NULL){
		ag->label=alab;
		ag->next = NULL;
		ag->prev = NULL;
		ag->exec = NULL;
		ag->pass = NULL;
		ag->S = NULL;
		ag->spp = NULL;
		ag->status = B_UNBOUND;
		ag->idx = agct++;ag->nbind=0;ag->ect=0;
		ag->biomass=0;
		ag->x=-1;
		ag->y=-1;
		return ag;
	}
	else{
		printf("mymalloc error\n");fflush(stdout);
		getchar();
		return NULL;
	}
}





/*******************************************************************************
* @brief print the agents
*
* @param[in] fp file pointer (including stdout)
*
* @param[in] spec either "NOW" or "NEXT" - nowhead or nexthead
*
* @param[in] verbose verbose output
*
* @return an agent, or NULL if not available
******************************************************************************/
void stringPM::AgentsPrint(FILE *fp, const char *spec, int verbose){

	s_ag *pag;
	pag = NULL;

	if(!strncmp("NOW",spec,strlen("NOW"))){
		pag = nowhead;
	}

	if(!strncmp("NEXT",spec,strlen("NEXT"))){
		pag = nexthead;
	}

	while(pag!=NULL){
		if(verbose){
			switch(pag->status){
				case B_UNBOUND:
					fprintf(fp,"Agent %6d,\texec=%4d\tnbind=%3d\tUNBOUND, %s\n",
							pag->idx,pag->ect,pag->nbind,pag->S);
					break;
				case B_ACTIVE:
					ReactionPrintState(fp,pag,pag->pass);
					break;
				default:
					break;
				}
		}
		else{
			switch(pag->status){
				case B_UNBOUND:
					fprintf(fp,"Agent %6d,\texec=%4d\tnbind=%3d\tUNBOUND, %s\n",
							pag->idx,pag->ect,pag->nbind,pag->S);
					break;
				case B_ACTIVE:
					fprintf(fp,"Agent %6d,\texec=%4d\tnbind=%3d\t ACTIVE, %s\n",
							pag->idx,pag->ect,pag->nbind,pag->S);
					fprintf(fp,"Agent %6d,\texec=%4d\tnbind=%3d\tPASSIVE, %s\n",
							pag->pass->idx,pag->pass->ect,pag->pass->nbind,
							pag->pass->S);
					//ReactionPrintState(stdout,pag,pag->pass);
					break;
				default:
					break;
				}
		}
		pag=pag->next;
	}
}




/******************************************************************************
 * @brief Print the state of agent with index 'idx'
 *
 * @details
 *
 * @param[in] fp file pointer, can be stdout for print to screen etc. NB no error checking for file stats
 *
 * @param[in] detail flags the level of detail
 *
 * @param[in] idx the index of the agent
 *
 * @return 0 regardless of succes (todo: fix this)
 *****************************************************************************/
int stringPM::AgentPrintWithIndex(FILE *fp, int detail, int idx){
	s_ag *pag;
	pag = nowhead;
	while(pag!=NULL){
		if(pag->idx == idx)
			switch(pag->status){
			case B_UNBOUND:
				printf("Agent %6d,\texec=%4d\tnbind=%3d\tUNBOUND, %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
				return 1;
			case B_ACTIVE:
				if(detail)
					ReactionPrintState(stdout,pag,pag->pass);
				else{
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tACTIVE,  %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tPASSIVE, %s\n",pag->pass->idx,pag->pass->ect,pag->pass->nbind,pag->pass->S);
				}
				//ReactionPrintState(stdout,pag,pag->pass);
				return 1;
			case B_PASSIVE:
				if(detail){
					ReactionPrintState(stdout,pag->exec,pag);
				}
				else{
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tACTIVE,  %s\n",pag->exec->idx,pag->exec->ect,pag->exec->nbind,pag->exec->S);
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tPASSIVE, %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
				}
				//ReactionPrintState(stdout,pag,pag->pass);
				return 1;
			}

		pag=pag->next;
	}
	return 0;
}





/******************************************************************************
 * @brief decay agent with a probability
 *
 * @param[in] pag the agent
 *
 * @return 1 if decay happens, 0 if not
 *****************************************************************************/
int stringPM::AgentAttemptDecay(s_ag *pag){

#ifdef LONG_DECAY
	//VARIABLE DECAY RATE BASED ON LENGTH OF STRING (ECAL 2009)
	float len = strlen(pag->S);
	float prob = 1./pow(len,2);//4./3.);
#else
	//CONSTANT DECAY RATE to match ECAL (ALife 2010 and on)
 	float prob = decayrate;//1./pow(65,2);//4./3.); //This is now done in load_decay...
#endif


	float rno = RandomBetween0And1();

#ifdef UNB_DECAY_ONLY
	if(rno<prob && pag->status == B_UNBOUND){
#else
	if(rno<prob && dodecay){
#endif
		//AgentUnbind(pag);
		AgentFree(pag);
		//TODO: sort out NULL-ing of freed agents
		//pag = NULL;
		return 1;
	}
	else
		return 0;
}







float stringPM::ReactionCalculateBindProbability(align *sw){
	float bprob = 0.;
	//This is the old bind prob, with a modifier for short strings:

	int l = sw->e1-sw->s1 < sw->e2-sw->s2 ? sw->e1-sw->s1 : sw->e2-sw->s2;
	if(l<=2)
		bprob=0;
	else{
		//bprob = pow(sw->score,l)/pow(l,l);
		//BRUTAL HACK:
		float s = sw->score<l-1.124? sw->score : l-1.124;
		bprob = s/(l-1.124);
	}

	return bprob;
}





/******************************************************************************
 * @brief align two molecular sequences using Smith-Waterman
 *
 * @param[in] a1 first agent
 *
 * @param[in] a2 second agent
 *
 * @param[in] sw struct to hold the alignment data
 *
 * @return the bind probability
 *****************************************************************************/
float stringPM::AgentsAlign(s_ag *a1, s_ag *a2, align *sw){

	float bprob;
	s_sw *swa;

	//SUGGEST: pass in pointer to the species - not its index
	swa = ReactionReadAlignmentFromSWList(swlist,a1->spp->spp,a2->spp->spp);

	if(swa==NULL){

		char *comp;
	
		comp = StringComplement(a1->S);

		/*bprob =*/ SmithWatermanAlignment(comp,a2->S,sw,blosum,0);

		free(comp);

		align sw2;

		/*bprob =*/ SmithWatermanAlignment(a1->S,a2->S,&sw2,blosum,0);

		//TODO: SUGGEST: pass in pointer to the species - not its index
		ReactionStoreAlignmentToSWList(&swlist,sw,a1->spp->spp,a2->spp->spp);
	}
	else{
		SmithWatermanDataFromAlignmentObject(swa,sw);
	}

	bprob = ReactionCalculateBindProbability(sw);


	if(verbose_bind){
		printf("Alignment:\nm1: %d to %d\nm2: %d to %d\nscore = %f\nProb = %f = %E\n",sw->s1,sw->e1,sw->s2,sw->e2,sw->score,bprob,bprob);
	}


	return bprob;
}





/******************************************************************************
 * @brief print spaces to line pointer up with position on string
 *
 * @param[in] fp file pointer
 *
 * @param[in] S the opcode string
 *
 * @param[in] p the pointer position
 *
 * @param[in] F toggle state of pointer - uppercase active, lowercase passive
 *
 * @param[in] c the pointer type (i,r,w,f)
 *****************************************************************************/
void stringPM::PointerPrintOffset(FILE *fp,const char *S,const char *p,int F, char c){
	int i,n=p-S;
	if(n<0){
		printf("Problem calculating pointer location\n");
		fflush(stdout);
	}
	for(i=0;i<n;i++)
		fprintf(fp," ");
	fprintf(fp,"%c\n",F?c-32:c);
}





/******************************************************************************
 * @brief print the current reaction state
 *
 * @param[in] fp file pointer
 *
 * @param[in] act active molecule
 *
 * @param[in] pas passive molecule
 *****************************************************************************/
void stringPM::ReactionPrintState(FILE *fp, s_ag *act, s_ag *pas){

	//Diagnostics to screen for passive:
	if(!strlen(pas->S))
		printf("Zero length passive string\n");
	fprintf(fp,"%6d:\n%s\n",pas->idx,pas->S);
	PointerPrintOffset(fp,pas->S,act->i[0],1-act->it,'i');
	PointerPrintOffset(fp,pas->S,act->f[0],1-act->ft,'f');
	PointerPrintOffset(fp,pas->S,act->r[0],1-act->rt,'r');
	PointerPrintOffset(fp,pas->S,act->w[0],1-act->wt,'w');

	//Diagnostics to screen for active:
	if(!strlen(act->S))
		printf("Zero length active string\n");
	fprintf(fp,"%6d:\n%s\n",act->idx,act->S);
	PointerPrintOffset(fp,act->S,act->i[1],act->it,'i');
	PointerPrintOffset(fp,act->S,act->f[1],act->ft,'f');
	PointerPrintOffset(fp,act->S,act->r[1],act->rt,'r');
	PointerPrintOffset(fp,act->S,act->w[1],act->wt,'w');

	act->len = strlen(act->S);
	pas->len = strlen(pas->S);


	if(pas->len<=(int) maxl){
		//printf("Passive string length = %d\n",pas->len);
	}
	else
		printf("Passive string length = %d - TOO LONG\n",pas->len);

	if(act->len<=(int) maxl){
		//printf("Active  string length = %d\n",act->len);
	}
	else
		printf("Active  string length = %d - TOO LONG\n",act->len);


}





/******************************************************************************
 * @brief Initialise the partners in a reaction

 * @param[in] A one agent
 *
 * @param[in] B the other agent
 *
 * @param[in] sw the SmithWaterman alignment data
 *
 * @return 1 if bind happens, 0 if not
 *****************************************************************************/
void stringPM::ReactionSetupExecution(s_ag *A, s_ag *B, align *sw){

	s_ag *active,*passive;
	int active_idx,passive_idx;

	if(sw->s1>=sw->s2){
		active = A;
		passive = B;
		active_idx = sw->s1;
		passive_idx = sw->s2;
	}
	else{
		active = B;
		passive = A;
		active_idx = sw->s2;
		passive_idx = sw->s1;
	}

	active->status = B_ACTIVE;
	active->exec = NULL;
	active->pass = passive;

	active->biomass = 0;

	passive->status = B_PASSIVE;
	passive->exec = active;
	passive->pass = NULL;

	active->f[0] = active->i[0] = active->r[0] = active->w[0] = &(passive->S[passive_idx]);//&(passive->S[sw->s1]);
	active->f[1] = active->i[1] = active->r[1] = active->w[1] = &(active->S[active_idx]);//&(active->S[sw->s2]);
	active->ft   = active->it   = active->rt   = active->wt = 1;

	passive->f[0] = passive->i[0] = passive->r[0] = passive->w[0] = 0;
	passive->f[1] = passive->i[1] = passive->r[1] = passive->w[1] = 0;
	passive->ft   = passive->it   = passive->rt   = passive->wt = 0;

#ifdef V_VERBOSE
	printf("Bind finished - looks like:\n");
	ReactionPrintState(stdout,active,passive);
#endif
}





/******************************************************************************
 * @brief Attempt to bind two molecules and initiate a reaction
 *
 * @details one agent has already been selected. Another is randomly chosen
 *          and a Smith-Waterman alignment is used to determine the chance
 *          of complementary binding
 *          todo: propensity is used, but I'm not sure if it always returns 1
 *
 * @param[in] pag
 *
 * @return 1 if bind happens, 0 if not
 *****************************************************************************/
int stringPM::AgentUnbind(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp){

	int found;
	int mass=0;

	if(pag->status==B_ACTIVE){
		mass = pag->biomass;
		pag->biomass = 0;
	}

	pag->status = B_UNBOUND;
	pag->pass = NULL;
	pag->exec = NULL;

	pag->ect=0;

	pag->f[0] = pag->i[0] = pag->r[0] = pag->w[0] = 0;
	pag->f[1] = pag->i[1] = pag->r[1] = pag->w[1] = 0;
	pag->ft   = pag->it   = pag->rt   = pag->wt = 0;

	found = SpeciesListUpdate(pag,sptype,update,pa,pp,mass);
	return found;
}





/******************************************************************************
 * @brief Attempt to bind two molecules and initiate a reaction
 *
 * @details one agent has already been selected. Another is randomly chosen
 *          and a Smith-Waterman alignment is used to determine the chance
 *          of complementary binding
 *          todo: propensity is used, but I'm not sure if it always returns 1
 *
 * @param[in] pag
 *
 * @return 1 if bind happens, 0 if not
 *****************************************************************************/
int stringPM::ReactionAttemptBind(s_ag *pag){

	int found=0;
	int count=0;
	float rno;
	s_ag *bag;
	align sw;
	float bprob;

	count = AgentsCount(nowhead,B_UNBOUND);

	//Calc propensity
	if(PropensityEquation(count)){
		found = 1;
	}

	if(found){
		bag = AgentSelectRandomly(nowhead,B_UNBOUND);
#ifndef BIND_ALL
		bprob = AgentsAlign(pag,bag,&sw);
#else
		bprob =1.0;
		sw.match = 1;		// the number of matching characters.
		sw.score = 1; 	// the score of the match
		sw.prob = 1.0;		// the probability of the match - used for determining events based on the score/match
		sw.s1=0;			// start of the match in string 1
		sw.e1=1;			// end of the match in string 1
		sw.s2=0;			// start of the match in string 2
		sw.e2=1;			// end of the match in string 2
#endif
		rno = RandomBetween0And1();
		if(rno<bprob){//Binding success!
			//figure out which is the executing string:
			ReactionSetupExecution(pag,bag,&sw);
			pag->nbind++;
			bag->nbind++;

			AgentAppend(&nexthead,pag);

			//uptime the second agent;
			AgentExtract(&nowhead,bag);
			AgentAppend(&nexthead,bag);
			energy--;
		}
		else{
			found = 0;
		}
	}
	PropensityRecord(count,found);
	return found;
}




/******************************************************************************
 * @brief pointer position relative to start of string
 *
 * @param[in] pag the chemical agent
 *
 * @param[in] headtype the pointer type
 *
 * @return the pointer position
 *****************************************************************************/
int stringPM::PointerPosition(s_ag *pag, char headtype){

	char *ph;
	char *ps;
	ps = NULL;
	ph = NULL;

	if(pag->status != B_ACTIVE)
		printf("ERROR: attempting headtype position for inactive string");

	switch(headtype){
	case 'w':
		ph = pag->w[pag->wt];
		if(pag->wt)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	case 'f':
		ph = pag->f[pag->ft];
		if(pag->ft)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	case 'i':
		ph = pag->i[pag->it];
		if(pag->it)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	case 'r':
		ph = pag->r[pag->rt];
		if(pag->rt)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	}

	return ph-ps;

}





/******************************************************************************
 * @brief reposition pointers if they are beyond the end of string post cleave
 *
 * @param[in] act the agent
 *
 * @return 0 always
 *****************************************************************************/
int stringPM::AgentRewindDanglingPtrs(s_ag* act){

	int plen,alen,pdist;
	char *ps;

	//PUT DANGLING POINTERS AT THE *END* OF THE STRINGS:
	//DO THE PASSIVE POINTERS FIRST:
	plen = strlen(act->pass->S);
	if(plen){
		ps = act->pass->S;

		pdist = act->i[0]-ps;
		if(pdist>plen || pdist<0)
			act->i[0]=ps+plen;

		pdist = act->r[0]-ps;
		if(pdist>plen || pdist<0)
			act->r[0]=ps+plen;

		pdist = act->w[0]-ps;
		if(pdist>plen || pdist<0)
			act->w[0]=ps+plen;

		pdist = act->f[0]-ps;
		if(pdist>plen || pdist<0)
			act->f[0]=ps+plen;
	}
	else{//Toggle everything off this string...
		act->i[0]=act->pass->S;
		act->r[0]=act->pass->S;
		act->w[0]=act->pass->S;
		act->f[0]=act->pass->S;

		act->it=1;
		act->rt=1;
		act->wt=1;
		act->ft=1;
	}


	//DO THE ACTIVE POINTERS NOW
	alen = strlen(act->S);
	ps = act->S;

	if(alen){
		pdist = act->i[1]-ps;
		if(pdist>alen || pdist<0)
			act->i[1]=ps+alen;

		pdist = act->r[1]-ps;
		if(pdist>alen || pdist<0)
			act->r[1]=ps+alen;

		pdist = act->w[1]-ps;
		if(pdist>alen || pdist<0)
			act->w[1]=ps+alen;

		pdist = act->f[1]-ps;
		if(pdist>alen || pdist<0)
			act->f[1]=ps+alen;
	}
	else{//Toggle everything off this string...
		act->i[1]=act->S;
		act->r[1]=act->S;
		act->w[1]=act->S;
		act->f[1]=act->S;

		if(plen){
			act->it=0;
			act->rt=0;
			act->wt=0;
			act->ft=0;
		}
	}

	//TODO: error checking on this!
	return 0;
}





/******************************************************************************
 * @brief handle zero-length strings (after a Cleave for example)
 *
 * @param[in] act the agent
 *
 * @return 0 if no problem; 1 if zero-length active string;
 *         2 zero-length passive string
 *****************************************************************************/
int stringPM::AgentCheckZeroLengthString(s_ag* act){

#ifdef VERBOSE
	ReactionPrintState(stdout,act,act->pass);
#endif
	//int len,pdist;
	//char *ps;

	//Sort pointers out first - even if there's going to be an error!
	AgentRewindDanglingPtrs(act);

	//Step 1: make sure act and pass *have* strings...
	if(!strlen(act->S)){
#ifdef VERBOSE
		printf("Zero length active string - dissoc\n");
#endif
		return 1;
		//if(!strlen(act->pass->S)){
		//	printf("Zero length active and passive strings - destroy\n");
		//	return 3;
		//}
	}
	if(!strlen(act->pass->S)){
#ifdef VERBOSE
		printf("Zero length passive string - dissoc\n");
#endif
		return 2;
	}



#ifdef VERBOSE
	ReactionPrintState(stdout,act,act->pass);
#endif

	return 0;

}



/******************************************************************************
 * @brief Copy operator "\="
 *
 * @details writes the symbol at the read pointer to the position of the
 *          write pointer, with mutation
 *
 * @param[in] act pointer to the active string (from which the partner string
 *            can be accessed)
 *
 * @return 0 if successful;
 *         -1 if attempt to write beyond maxl;
 *         -2 if attempt to read beyond maxl;
 *****************************************************************************/
int stringPM::OpcodeCopy(s_ag *act){

	//s_ag *pass;
	//pass = act->pass;
	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..


	//MUTATION RATES:
	//THESE ARE HARD-CODED FOR NOW - THEY SHOULD BE DERIVED FROM THE BLOSUM SOMEHOW...
	//const float indelrate = 0.0005,		subrate=0.375;//0.0749/2
	//const float indelrate = 0.0000306125,	subrate=0.01;//0.02
	//const float indelrate = 0.00006125,	subrate=0.05;//0.02
	//const float indelrate = 0.000125,		subrate=0.1;//0.02
	//const float indelrate = 0.000005,		subrate=0.0375;//0.0749/2
	//const float indelrate = 0., 			subrate = 0.;

	if(!domut){
		indelrate = subrate =0;
	}

	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	
	if( (PointerPosition(act,'w'))>=(int) maxl){

//#ifdef DODEBUG
//		printf("Write head out of bounds: %d\n",p);
//#endif

		//just to make sure no damage is done:
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';

		act->i[act->it]++;

		return -1;
	}
	if(PointerPosition(act,'r')>=(int) maxl){
		printf("Read head out of bounds\n");

		act->i[act->it]++;

		return -2;
	}
	//TODO: see note on error check in "another" hcopy below

	if(*(act->r[act->rt]) == 0){
		//possibly return a negative value and initiate a b
		safe = 0;
		//return -3;
	}
	//if(h_pos(act,'w')>=maxl){
	//	//We are at the end of a copy...
	//	//...so just increment *R
	//	act->r[act->rt]++;
	//	safe = 0;
	//}

	if(safe){
		rno=RandomBetween0And1();
		if(rno<indelrate){//INDEL

			//should follow the blosum table for this....
			rno=RandomBetween0And1();
			if(rno<0.5){//insert
				//first do a straight copy..
				*(act->w[act->wt])=*(act->r[act->rt]);

				//no need to test for granular here since we are inserting...
				act->w[act->wt]++;

				//Then pick a random instruction:
				cidx = (float) RandomBetween0And1() * blosum->N;

				//insert the random instruction
				*(act->w[act->wt])=blosum->key[cidx];
				if(granular_1==0){
					act->w[act->wt]++;
				}
			}
			else{//delete
				act->i[act->it]++;
			}

			if(granular_1==0){
				act->r[act->rt]++;
			}

		}
		else{
			if(rno<subrate+indelrate){//INCREMENTAL MUTATION
				cidx = OpcodeAdjacent(*(act->r[act->rt]),blosum);
				*(act->w[act->wt])=cidx;
			}
			else{//NO MUTATION
				*(act->w[act->wt])=*(act->r[act->rt]);
			}
			if(granular_1==0){
				act->w[act->wt]++;
				act->r[act->rt]++;
			}
		}
	}
	//update lengths
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);


	act->i[act->it]++;


#ifdef VERBOSE
	if(mut)
	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
	act->biomass++;
	biomass++;
	return 0;
}



/******************************************************************************
 * @brief cleave
 *
 * @param[in] act
 *
 * @return which agents have been destroyed (if any) 0: none; 1: active only
 *         2: passive only; 3: both
 *****************************************************************************/
int stringPM::OpcodeCleave(s_ag *act){

	int dac = 0,cpy;
	s_ag *c,*pass,*csite;

	pass = act->pass;

	//pick the mol containing the cleave site:
	csite = act->ft?act:pass;

	if(act->f[act->ft]-csite->S < csite->len){

		//1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT

		//Can't really say what the label is easily - for ECAL, it's always pass
		c = AgentMake(pass->label);//,1);

		//Copy the cleaved string to the agent
		char *cs;

		c->S =(char *) malloc(maxl0*sizeof(char));
		memset(c->S,0,maxl0*sizeof(char));

		cs = csite->S;
		cpy = strlen(cs);
		//Check that we aren't creating a zero-length molecule:
		if(!cpy){
			printf("WARNING: Zero length molecule being created!\n");
		}

		//Make the parent structure: ALL DONE NOW IN SpeciesListUpdate
		//c->pp = splist->ParentsMake(act->spp,pass->spp);

		cpy -= act->f[act->ft]-cs;

		strncpy(c->S,act->f[act->ft],cpy);
		c->len = strlen(c->S);
#ifdef VERBOSE
		printf("String %d created:\n%s\n",c->idx,c->S);
#endif

		//ALL DONE NOW IN SpeciesListUpdate
		//Fill in the birth certificate:
		//Parents now set in SpeciesListUpdate
		//c->paspp = act->spp;
		//c->ppspp = pass->spp;

		//Check the lineage
		SpeciesListUpdate(c,'C',1,act->spp,pass->spp,act->biomass);
		act->biomass=0; //reset this; we might continue to make stuff!

		//append the agent to nexthead
		AgentAppend(&nexthead,c);

		//2: HEAL THE PARENT

		memset(act->f[act->ft],0,cpy*sizeof(char));

		csite->len = strlen(csite->S);

		if((dac = AgentCheckZeroLengthString(act))){
			switch(dac){
			case 1://Destroy active - only append passive
				AgentUnbind(pass,'P',1,act->spp,pass->spp);
				AgentAppend(&nexthead,pass);
				AgentFree(act);
				act = NULL;
				break;
			case 2://Destroy passive - only append active
				AgentUnbind(act,'A',1,act->spp,pass->spp);
				AgentAppend(&nexthead,act);
				AgentFree(pass);
				pass = NULL;
				break;
			case 3://Destroy both
				printf("This should never happen\n");
				AgentUnbind(act,'A',1,act->spp,pass->spp);
				AgentUnbind(pass,'P',1,act->spp,pass->spp);
				AgentFree(act);
				act = NULL;
				AgentFree(pass);
				pass = NULL;
				break;
			default://This can't be right can it? - NB dac = 0 covered here - make explicit!
				if(act->ft == act->it){
					act->i[act->it]--;
				}
				break;
			}
		}
	}
	if(!dac){
		act->i[act->it]++;
		//dac=-1;
	}

	return dac;
}



/******************************************************************************
 * @brief execute the current opcode in a reaction
 *
 * @param[in] act the active string
 *
 * @param[in] pass the passive string
 *
 * @return 1 always - to indicate the reaction has changed
 *         todo: maybe some error handling here would be good!
 *****************************************************************************/
int stringPM::ReactionExecuteOpcode(s_ag *act, s_ag *pass){

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
		tmp = OpcodeSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
		act->f[act->ft] = tmp;
		act->i[act->it]++;
		break;

		/*
		//put the active flow head on the string where the active ip is:
		act->ft = act->it;
		char *cs;
		if(act->it)
			cs = act->S;
		else
			cs = act->pass->S;
		act->f[act->ft] = HSearch(act->i[act->it],cs,blosum);
		act->i[act->it]++;
		break;
		*/

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
		if(OpcodeCopy(act)<0){
			AgentUnbind(act,'A',1,act->spp,pass->spp);
			AgentUnbind(pass,'P',1,act->spp,pass->spp);
		}
		break;


	/************
	 *   INC_R  *
	 ************/
	case '+'://h-copy
		if(granular_1==1){
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
			act->i[act->it]=OpcodeIf(act->i[act->it],act->r[act->rt],act->S,blosum,maxl);
			break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			if((/*dac =*/ OpcodeCleave(act) )){

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
			AgentUnbind(act,'A',1,act->spp,pass->spp);
			AgentUnbind(pass,'P',1,act->spp,pass->spp);

			break;

	default://Just increment the i-pointer
		act->i[act->it]++;
		break;
	}
#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	ReactionPrintState(stdout,act,pass);
#endif


	if(safe_append){
		act->ect++;
		AgentAppend(&nexthead,act);
		AgentAppend(&nexthead,pass);
	}
	energy--;


	return 1;

}





/******************************************************************************
 * @brief process all agents in this timestep
 *
 * @detail While there's some agents left in "now"
 *      -pick an agent at random
 *      -if there are any rules applying to the agent
 *          -find the maximum number of rules on the lhs of all rules applying to the agent = maxl
 *          switch(maxl)
 *          -case 2: We have a bonding rule to test
 *          	-go through all other agents.
 *          		-If closer than the current closest agent && there's a bonding rule for the two
 *						-set 2nd agent to agent b.
 *						-update closest
 *				- if we have a reactive pair
 *					-if we've enough energy and a random number is less than reaction prob
 *						-do the reaction and put the product in the 'next' bucket
 * 			-case 1: Dissocn or decay rule to test (assuming only one of these rules)
 *				- get the appropriate rule and find if its dissoc or decay
 *				- if dissoc
 *					- do the dissoc
 *				- else
 *					- do the decay
 *		- if no reaction has been carried out
 *			- move the agent to the next bucket
 *****************************************************************************/
void stringPM::TimestepIncrement(){
	s_ag *pag;

	//SUGGEST: write function to count what's around (saves rechecking every time)
	//countstates();

	//SUGGEST: IF WE ARE PUTTING ENERGY IN FROM OUTSIDE, AND NO ENERGY CAN BE PRODUCED
	//while(nowhead!=NULL && energy){
	//}
	//if(nowhead != NULL){
	//	pag = nowhead;
	//	append_ag(pag); //check that pag->next is also appended
	// NB: will have to have a different scheme for decay, since that doesn't require energy - sample from a binomial?

	while(nowhead!=NULL){

		pag = AgentSelectRandomly(nowhead,-1);
		AgentExtract(&nowhead,pag);

		//TODO: This is a debug option:
		//if(0){//extit>=10){
		//	if(pag->status==B_PASSIVE)
		//		//if(pag->exec->idx==1041){
		//			ReactionPrintState(stdout,pag->exec,pag);
		//		//}
		//	if(pag->status==B_ACTIVE){
		//		//if(pag->idx==1041)
		//			ReactionPrintState(stdout,pag,pag->pass);
		//	}
		//	fflush(stdout);
		//}

		//extract any partner:
		s_ag *bag;
		bag = NULL;

		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			//if(pag->idx == 214 && extit > 79)
			//	printf("Active 214\n");
			bag = pag->pass;
			AgentExtract(&nowhead,bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			AgentExtract(&nowhead,bag);
			break;
		}

		//TODO: This looks like a debugging step - delete?
		if(pag->label=='P' && pag->status != B_UNBOUND  ){
			if(AgentPrintWithIndex(stdout,1,pag->idx))
				fflush(stdout);
		}

		int dc = AgentAttemptDecay(pag);
		if(dc){//we must check what else needs to be destroyed...
			if(bag!=NULL){
				AgentFree(bag);
				bag = NULL;
			}
		}
		else{
			int changed = 0;
			if(energy>0){
				switch(pag->status){
				case B_UNBOUND:
					//seek binding partner, set binding states.
					changed = ReactionAttemptBind(pag);
					//if(!changed)
					//	changed = testdecay(pag);

					break;
				case B_PASSIVE:

					//extract_ag(&nowhead,pag->exec);
					changed = ReactionExecuteOpcode(pag->exec,pag);

					break;
				case B_ACTIVE:

					//extract_ag(&nowhead,pag->pass);
					changed = ReactionExecuteOpcode(pag,pag->pass);
					break;
				default:
					printf("ERROR: agent with unknown state encountered!\n");
				}
			}
			if(!changed){
				AgentAppend(&nexthead,pag);
				if(bag!=NULL)
					AgentAppend(&nexthead,bag);

			}
		}
	}

	//This used to be update()...
	nowhead = nexthead;
	nexthead = NULL;
}





/*
//TODO: see where this differs from other versions of hcopy! 
int stringPM::speig_hcopy(s_ag *act){

	//s_ag *pass;
	//pass = act->pass;
	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..

	if(!domut){
		indelrate = subrate =0;
	}

	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	int p;
	if( (p = h_pos(act,'w'))>=(int) maxl){
		printf("Write head out of bounds: %d\n",p);
		//just to make sure no damage is done:
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';

		act->i[act->it]++;
		safe = 0;
		return -1;
	}
	if(h_pos(act,'r')>=(int) maxl){
		printf("Read head out of bounds\n");
		act->i[act->it]++;
		safe = 0;
		return -2;
	}




	if(*(act->r[act->rt]) == 0){
		//possibly return a negative value and initiate a b
		safe = 0;
		//return -3;
	}

	if(safe){

		//see if we are overwriting or not:
		int rm,wm=-1;
		if(*(act->w[act->wt])){
			wm=tab_idx(*(act->w[act->wt]),blosum);
		}

		const float speig_idrate = 0.001;
		float winc=rand0to1();
		float rinc=rand0to1();

		//todo: make sure no increments happen if the symbol (or mutant) is not available

		rno=rand0to1();
		if(rno<subrate){//INCREMENTAL MUTATION

			cidx = sym_from_adj(*(act->r[act->rt]),blosum);
			rm = tab_idx(cidx,blosum);
			if(mass[rm]){
				*(act->w[act->wt])=cidx;
				if(winc>speig_idrate)
					act->w[act->wt]++;

				if(rinc>speig_idrate)
					act->r[act->rt]++;//possible deletion here...

				if(!(wm<0)){
					mass[wm]++;
				}
				mass[rm]--;
			}
		}
		else{//NO MUTATION (but possible sub via comass effects)
			//cidx = sym_from_adj(*(act->r[act->rt]),blosum);
			rm = tab_idx(*(act->r[act->rt]),blosum);
			if(mass[rm]){
				*(act->w[act->wt])=*(act->r[act->rt]);

				if(winc>speig_idrate)
					act->w[act->wt]++;

				if(rinc>speig_idrate)
					act->r[act->rt]++;

				if(!(wm<0)){
					mass[wm]++;
				}
				mass[rm]--;
			}
			else{
				cidx = sym_from_adj(*(act->r[act->rt]),blosum);
				rm = tab_idx(cidx,blosum);
				if(mass[rm]){
					if(winc>speig_idrate)
						act->w[act->wt]++;

					if(rinc>speig_idrate)
						act->r[act->rt]++;

					act->w[act->wt]++;
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
			}
		}

	}
	//update lengths
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);
	act->i[act->it]++;

#ifdef VERBOSE
	if(mut)
	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
	act->biomass++;
	biomass++;
	return 0;
}*/


/////////////////////////////////////////////////////////start of comass stuff



int stringPM::comass_hcopy(s_ag *act){

	//s_ag *pass;
	//pass = act->pass;
	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..


	//MUTATION RATES:
	//THESE ARE HARD-CODED FOR NOW - THEY SHOULD BE DERIVED FROM THE BLOSUM SOMEHOW...
	//const float indelrate = 0.0005,subrate=0.375;//0.0749/2
	//const float indelrate = 0.0000306125,subrate=0.01;//0.02
	//const float indelrate = 0.00006125,subrate=0.05;//0.02
	//const float indelrate = 0.000125,subrate=0.1;//0.02
	//const float indelrate = 0.000005,subrate=0.0375;//0.0749/2

	//const float indelrate = 0., 		subrate = 0.;

	if(!domut){
		indelrate = subrate =0;
	}
	//Make sure the recorded lengths are current
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	int p;

	//Check positions of the pointers
	//TODO: shouldn't be casting maxl to int - but can h_pos ever return negative?
	if( (p = PointerPosition(act,'w'))>=(int) maxl){
		printf("Write head out of bounds: %d\n",p);
		//just to make sure no damage is done:
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';

		act->i[act->it]++;
		return -1;
	}
	if(PointerPosition(act,'r')>=(int) maxl){
		printf("Read head out of bounds\n");
		act->i[act->it]++;
		return -2;
	}
	//TODO: Are we handling the above errors ok?
	//Check that we aren't off the end of the string, but within the allocated memory:
	if(*(act->r[act->rt]) == 0){
		//possibly return a negative value and initiate a b
		safe = 0;
		//return -3;
	}


	if(safe){

		//see if we are overwriting or not:
		int rm,wm=-1;
		if(*(act->w[act->wt])){
			wm=OpcodeIndex(*(act->w[act->wt]),blosum);
		}


		rno=RandomBetween0And1();
		if(rno<indelrate){//INDEL - we should never be doing this in comass!

			//should follow the blosum table for this....
			rno=RandomBetween0And1();
			if(rno<0.5){//insert

				//first do a straight copy..
				*(act->w[act->wt])=*(act->r[act->rt]);
				act->w[act->wt]++;

				//Then pick a random instruction:
				rm = (float) RandomBetween0And1() * blosum->N;

				//Check there's mass for this symbol:
				if(mass[rm]){
					//insert the random instruction
					*(act->w[act->wt])=blosum->key[rm];
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
				act->w[act->wt]++;
			}
			else{//delete
				act->i[act->it]++;

				//simply increment the read head without doing anything else
			}
			act->r[act->rt]++;
			//act->i[act->it]++;
		}
		else{
			if(rno<subrate+indelrate){//INCREMENTAL MUTATION - we should never be doing this either

				cidx = OpcodeAdjacent(*(act->r[act->rt]),blosum);
				rm = OpcodeIndex(cidx,blosum);
				if(mass[rm]){
					*(act->w[act->wt])=cidx;
					act->w[act->wt]++;
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
				else{//
					//simply increment the read head without doing anything else
				}
				act->r[act->rt]++;//possible deletion here...
			}
			else{//NO MUTATION (but possible sub via comass effects)

				rm = OpcodeIndex(*(act->r[act->rt]),blosum);
				if(mass[rm]){
					*(act->w[act->wt])=*(act->r[act->rt]);
					act->w[act->wt]++;
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
				else{
					cidx = OpcodeAdjacent(*(act->r[act->rt]),blosum);
					rm = OpcodeIndex(cidx,blosum);
					if(mass[rm]){
						*(act->w[act->wt])=cidx;
						act->w[act->wt]++;
						if(!(wm<0)){
							mass[wm]++;
						}
						mass[rm]--;
					}
				}
				act->r[act->rt]++;
			}
		}
	}
	//update lengths
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);
	act->i[act->it]++;

#ifdef VERBOSE
	if(mut)
	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
	act->biomass++;
	biomass++;
	return 0;
}



int stringPM::comass_ReactionExecuteOpcode(s_ag *act, s_ag *pass){

	int finished=0;
	char *tmp;
	int safe_append=1;

	switch(*(act->i[act->it])){

	case '$'://h-search

		char *cs;
		if(act->ft)
			cs = act->S;
		else
			cs = act->pass->S;
		tmp = OpcodeSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
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
		if(comass_hcopy(act)<0){
			AgentUnbind(act,'A',1,act->spp,pass->spp);
			AgentUnbind(pass,'P',1,act->spp,pass->spp);
			finished = 1;
		}
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
			act->i[act->it]=OpcodeIf(act->i[act->it],act->r[act->rt],act->S,blosum,maxl);
			break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			if((/*dac =*/ OpcodeCleave(act) )){
				//AgentUnbind(act);
				//AgentUnbind(pass);
				finished = 1;
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
			AgentUnbind(act,'A',1,act->spp,pass->spp);
			AgentUnbind(pass,'P',1,act->spp,pass->spp);

			finished = 1;
			break;

	default://Just increment the i-pointer
		act->i[act->it]++;
		break;
	}
#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	ReactionPrintState(stdout,act,pass);
#endif


	if(safe_append){
		act->ect++;
		AgentAppend(&nexthead,act);
		AgentAppend(&nexthead,pass);
	}
	energy--;

	//TODO: propagate this up if needed, or get rid...
	return finished;
}



//This to be called AFTER agents and blosum have been loaded.
//Sets the mass for everything to a single value, read from config.
int stringPM::load_comass(const char *fn, int verbose){
	//int massval = 2000;
	FILE *fp;
	unsigned int massval;

	mass = (int*) malloc(blosum->N * sizeof(int));

	if((fp=fopen(fn,"r"))!=NULL){

		int finderr=ParameterReadUnsignedInt(fp,"MASS",&massval,verbose);

		report_param_error(finderr, 1);


		fclose(fp);
		//Load the massvalue into the mass table;
		for(int i=0;i<blosum->N;i++){
			mass[i]=massval;
		}
		s_ag *pag;
		for(pag=nowhead;pag!=NULL;pag=pag->next){
			update_mass(pag->S,strlen(pag->S),-1, 1);
		}
	}
	return 0;
}


int stringPM::set_mass(const int *param){

	int i;

	mass = (int*) malloc(blosum->N * sizeof(int));

	for(i=0;i<blosum->N;i++){
		mass[i]=param[i];
	}
	s_ag *pag;
	for(pag=nowhead;pag!=NULL;pag=pag->next){
		if(update_mass(pag->S,strlen(pag->S),-1, 1)){
			//SpeciesListUpdate(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp, int mass){
			SpeciesListUpdate(pag,'I',1,NULL,NULL,0);
		}
	}

	return 0;
}





int stringPM::update_mass(char *S, int len, int val, const int doconcat){
	
	int updated=0;
	for(int i=0;i<len;i++){
		int c = OpcodeIndex(S[i],blosum);
		mass[c] += val;
		//Strategy: if there isn't the mass available for the string, delete the code
		if(doconcat && mass[c]<0){
			mass[c] -= val;
			for(int j = i+1; j< len; j++){
				S[j-1] = S[j];
			}
			S[--len]=0;
			i--; //decrement i so that we are checking the left-shifted character next
			updated=1;
		}
	}
	return updated;
}


int stringPM::comass_free_ag(s_ag *pag){

	if(pag->S != NULL){
		update_mass(pag->S,strlen(pag->S),1,0);//TODO: the last argument specifies whether to concatenate if no mass available, but since we are adding, it'll never be called
		//printf("destroying agent %d, code = %s\n",pag->idx,pag->S);
		free(pag->S);
	}

	free(pag);
	//pag = NULL;

	return 0;
}


int stringPM::comass_AgentAttemptDecay(s_ag *pag){

#ifdef LONG_DECAY
	//VARIABLE DECAY RATE BASED ON LENGTH OF STRING (ECAL 2009)
	float len = strlen(pag->S);
	float prob = 1./pow(len,2);//4./3.);
#else
	//CONSTANT DECAY RATE to match ECAL (ALife 2010 and on)
	float prob = 1./pow(65,2);//4./3.);
#endif


	float rno = RandomBetween0And1();

#ifdef UNB_DECAY_ONLY
	if(rno<prob && pag->status == B_UNBOUND){
#else
	if(rno<prob && dodecay){
#endif
		//AgentUnbind(pag);
		comass_free_ag(pag);
		return 1;
	}
	else
		return 0;
}

void stringPM::comass_TimestepIncrement(){
	s_ag *pag;


	//SUGGEST: write function to count what's around (saves rechecking every time)
	//countstates();

	//SUGGEST: IF WE ARE PUTTING ENERGY IN FROM OUTSIDE, AND NO ENERGY CAN BE PRODUCED
	//while(nowhead!=NULL && energy){
	//}
	//if(nowhead != NULL){
	//	pag = nowhead;
	//	append_ag(pag); //check that pag->next is also appended
	// NB: will have to have a different scheme for decay, since that doesn't require energy - sample from a binomial?

	while(nowhead!=NULL){

		pag = AgentSelectRandomly(nowhead,-1);
		AgentExtract(&nowhead,pag);

		//if(0){//extit>=10){
		//	if(pag->status==B_PASSIVE)
		//		//if(pag->exec->idx==1041){
		//			ReactionPrintState(stdout,pag->exec,pag);
		//		//}
		//	if(pag->status==B_ACTIVE){
		//		//if(pag->idx==1041)
		//			ReactionPrintState(stdout,pag,pag->pass);
		//	}
		//	fflush(stdout);
		//}

		//extract any partner:
		s_ag *bag;
		bag = NULL;

		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			//if(pag->idx == 214 && extit > 79)
			//	printf("Active 214\n");
			bag = pag->pass;
			AgentExtract(&nowhead,bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			AgentExtract(&nowhead,bag);
			break;
		}

		if(pag->label=='P' && pag->status != B_UNBOUND  ){
			if(AgentPrintWithIndex(stdout,1,pag->idx))
				fflush(stdout);
		}

		int dc = comass_AgentAttemptDecay(pag);
		if(dc){//we must check what else needs to be destroyed...
			if(bag!=NULL){
				comass_free_ag(bag);
			}
		}
		else{
		
			int changed = 0;
			if(energy>0){
				switch(pag->status){
				case B_UNBOUND:
					//seek binding partner, set binding states.
					changed = ReactionAttemptBind(pag);
					//if(!changed)
					//	changed = testdecay(pag);

					break;
				case B_PASSIVE:

					//extract_ag(&nowhead,pag->exec);
					changed = comass_ReactionExecuteOpcode(pag->exec,pag);

					break;
				case B_ACTIVE:

					//extract_ag(&nowhead,pag->pass);
					changed = comass_ReactionExecuteOpcode(pag,pag->pass);
					break;
				default:
					printf("ERROR: agent with unknown state encountered!\n");
				}
			}
			if(!changed){
				AgentAppend(&nexthead,pag);
				if(bag!=NULL)
					AgentAppend(&nexthead,bag);

			}
		}
	}

	nowhead = nexthead;
	nexthead = NULL;
}





/////////////////////////////////////////////////////////end of comass stuff




/////////////////////////////////////////////////////////start of energetic stuff


int stringPM::energetic_exec_step(s_ag *act, s_ag *pass){//pset *p,char *s1, swt *T){

	int finished=0;
	char *tmp;
	int safe_append=1;

	if(energy>0){
			switch(*(act->i[act->it])){//*iptr[it]){

				case '$'://h-search
					//act->ft = act->it;
					char *cs;
					if(act->ft)
						cs = act->S;
					else
						cs = act->pass->S;
					tmp = OpcodeSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
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
					if(OpcodeCopy(act)<0){
						AgentUnbind(act,'A',1,act->spp,pass->spp);
						AgentUnbind(pass,'P',1,act->spp,pass->spp);
						finished = 1;
					}
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
						act->i[act->it]=OpcodeIf(act->i[act->it],act->r[act->rt],act->S,blosum,maxl);
						break;


				/************
				 *  CLEAVE  *
				 ************/
				case '%':
						if((/*dac =*/ OpcodeCleave(act) )){
							finished = 1;
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
						AgentUnbind(act,'A',1,act->spp,pass->spp);
						AgentUnbind(pass,'P',1,act->spp,pass->spp);

						finished = 1;
						break;

				default://Just increment the i-pointer
					act->i[act->it]++;
					break;
		}
	}
	else{//Just increment the i-pointer
		act->i[act->it]++;
	}

#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	ReactionPrintState(stdout,act,pass);
#endif


	if(safe_append){
		act->ect++;
		AgentAppend(&nexthead,act);
		AgentAppend(&nexthead,pass);
	}
	if(energy>0)
		energy--;



	//return 1;
	return finished;
}



int stringPM::energetic_attempt_bind(s_ag *pag){

	int found=0;
	int count=0;
	const int nrglim = 3000;
	float
		rno;
	s_ag *bag;
	align sw;
	float bprob;

	count = AgentsCount(nowhead,B_UNBOUND);

	//Calc propensity
	if(PropensityEquation(count)){
		found = 1;
	}

	if(found){
		bag = AgentSelectRandomly(nowhead,B_UNBOUND);
#ifndef BIND_ALL
		bprob = AgentsAlign(pag,bag,&sw);
#else
		bprob =1.0;
		sw.match = 1;		// the number of matching characters.
		sw.score = 1; 	// the score of the match
		sw.prob = 1.0;		// the probability of the match - used for determining events based on the score/match
		sw.s1=0;			// start of the match in string 1
		sw.e1=1;			// end of the match in string 1
		sw.s2=0;			// start of the match in string 2
		sw.e2=1;			// end of the match in string 2
#endif
		rno = RandomBetween0And1();
		if(rno<bprob && ((1-bprob)*nrglim < energy   ) ){//Binding success!
			//figure out which is the executing string:
			ReactionSetupExecution(pag,bag,&sw);
			pag->nbind++;
			bag->nbind++;

			AgentAppend(&nexthead,pag);

			//uptime the second agent;
			AgentExtract(&nowhead,bag);
			AgentAppend(&nexthead,bag);
			energy--;
		}
		else{
			found = 0;
		}
	}
	PropensityRecord(count,found);
	return found;
}


void stringPM::energetic_TimestepIncrement(){
	s_ag *pag;

	while(nowhead!=NULL){

		pag = AgentSelectRandomly(nowhead,-1);
		AgentExtract(&nowhead,pag);

		//extract any partner:
		s_ag *bag;
		bag = NULL;
		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			//if(pag->idx == 214 && extit > 79)
			//	printf("Active 214\n");
			bag = pag->pass;
			AgentExtract(&nowhead,bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			AgentExtract(&nowhead,bag);
			break;
		}

		if(pag->label=='P' && pag->status != B_UNBOUND  ){
			if(AgentPrintWithIndex(stdout,1,pag->idx))
				fflush(stdout);
		}

		int dc = AgentAttemptDecay(pag);
		if(dc){//we must check what else needs to be destroyed...
			if(bag!=NULL){
				AgentFree(bag);
				bag = NULL;
			}
		}
		else{
			int changed = 0;

			switch(pag->status){
			case B_UNBOUND:
				if(energy>0){
					//seek binding partner, set binding states.
					changed = energetic_attempt_bind(pag);
					//if(!changed)
					//	changed = testdecay(pag);
				}

				break;
			case B_PASSIVE:

				//extract_ag(&nowhead,pag->exec);
				changed = energetic_exec_step(pag->exec,pag);

				break;
			case B_ACTIVE:

				//extract_ag(&nowhead,pag->pass);
				changed = energetic_exec_step(pag,pag->pass);
				break;
			default:
				printf("ERROR: agent with unknown state encountered!\n");
			}


			if(!changed){
				AgentAppend(&nexthead,pag);
				if(bag!=NULL)
					AgentAppend(&nexthead,bag);

			}
		}
	}

	nowhead = nexthead;
	nexthead = NULL;
}
/////////////////////////////////////////////////////////end of energetic stuff












/*******************************************************************************
 * @brief move "now" to "next" by pointer assignment
 *
 * @details nowhead = nexthead; nexthead = NULL;
 ******************************************************************************/
void stringPM::UpdateNowNext(){

	nowhead = nexthead;
	nexthead = NULL;
}





/*******************************************************************************
 * @brief free the swt object
 *
 * @param[in] pSWT pointer to the swt object
 *
 * @param[in] verbose flag
 ******************************************************************************/
void stringPM::SmithWatermanFree(swt *pSWT, int verbose){


	//print for debug.
	int i;
	if(verbose){
		printf("   ");
		for(i=0;i<pSWT->N;i++)
			printf("%c   ",pSWT->key[i]);
		printf("\n");
		for(i=0;i<pSWT->N;i++){
			int j;
			
			printf("%c  ",pSWT->key[i]);
			for(j=0;j<pSWT->N;j++){
				printf("%0.2f ",pSWT->T[i][j]);
			}
			printf("\n");
		}
		//Print the indel value
		printf("-   ");
		for(i=0;i<pSWT->N;i++){
			printf("%02d ",(int) pSWT->T[pSWT->N][i]);
		}
		fflush(stdout);
	}

	//Free elements of pSWT
	if(pSWT->T != NULL){
	
		free(pSWT->key);
		pSWT->key= NULL;
		
		for(i=0;i<pSWT->N+1;i++){
			float *pT;
			pT = pSWT->T[i];
			free(pT);
		}
		free(pSWT->T);
		pSWT->T=NULL;

		for(i=0;i<pSWT->N;i++){
			int *pI;
			pI = pSWT->adj[i];
			free(pI);
		}
		free(pSWT->adj);
		pSWT->adj = NULL;
	}


}





/*******************************************************************************
 * @brief free the grid
 ******************************************************************************/
void stringPM::GridFree(){

	if(grid !=NULL){

		for(int i=0;i<grid->gridx;i++){
			free(grid->status[i]);
			free(grid->grid[i]);
		}
		free(grid->status);
		free(grid->grid);
		free(grid);
		grid = NULL;
	}
}



// TODO(sjh): this should (should this?) be the destructor! ~stringPM
/*******************************************************************************
* @brief clear the structures in a run so the run can be repeated
*
* @param[in] verbose verbose output
*******************************************************************************/
void stringPM::BucketReset(int verbose){

	s_ag	*agp,*agp2;

	if(verbose){
		printf("Starting spatial clearout..");fflush(stdout);
	}

	SmithWatermanFree(blosum,verbose);

	agp = nowhead;
	while(agp!=NULL){
		agp2=agp->next;
		AgentFree(agp);
		agp=agp2;
	}

	agp=nexthead;
	while(agp!=NULL){
		agp2=agp->next;
		AgentFree(agp);
		agp=agp2;
	}

	GridFree();

	SmithWatermanListFree(&swlist);

	SetHeadsAndDefaults();

}





//todo(sjh): written to check buffer overruns on strings - potentially useful!
/*
void stringPM::sanity_check(){
	s_ag *pag;
	pag = nowhead;
	while(pag!=NULL){

		if(pag->S[maxl]!=0){
			printf("Problem with string %d, status ",pag->idx);
			printf(" out of range: %s", pag->S);

		}

		pag=pag->next;

	}
}*/

//SPECIES ANALYSIS FUNCTIONS






/*******************************************************************************
* @brief update the list of species with an agent
*
* @details This is called from CLEAVE, otherwise we can't tell if its in
*          the middle of being constructed....
*
* @param[in] p the agent
*
* @param[in] sptype the 'type' of species
*
* @param[in] add flag to say whether to add to the list
*
* @param[in] paspp pointer to active species list (???)
*
* @param[in] ppspp pointer to passive species list (???)
*
* @param[in] mass the number of characters in the string (???)
*
* @return 0 regardless of succes (todo: fix this)
*******************************************************************************/
int stringPM::SpeciesListUpdate(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp, int mass){

	l_spp *sp;
	sp = spl->species_list;
	int found=0;
	int novel=0;

	while(sp!=NULL&&!found){
		if(!strcmp(sp->S,p->S)){//we've found a match on the string
			found=sp->spp;
			if(add){
				//Need to check whether this is a dissociating partner in a reaction that has changed:
				if(p->spp!=NULL){//We haven't decided what the spp is yet
					if(p->spp->spp != sp->spp){
						sp->count++;
						novel=1;
					}
					//novel=0 IF dissociating and no new spp are produced.
				}
				else{
					sp->count++;
					novel=1;
				}
			}
			break;
		}
		sp = sp->next;
	}

	if(add){//Only do this if we are adding to the list (not just checking reaction-space
		if(!found){
			sp = spl->SpeciesMakeFromAgent(p,timestep,maxl0);

			sp->sptype = sptype;
			sp->biomass += mass;
			spl->SpeciesPrependToList(sp); //append_lspp(sp);
			novel=1;
		}

		//Now sort out the parentage:
		p->spp = sp;//->spp;
		if(novel)
			p->pp = spl->ParentsFindOrMake(p->spp, paspp, ppspp);

	}

	return found;
}





/*******************************************************************************
* @brief update the list of species with an agent
*
* @param[in] n the index (id)
*
* @return the species, or NULL if not found;
*******************************************************************************/
l_spp * stringPM::SpeciesFromListByIndex(int n){
	l_spp * sp;
	sp = spl->species_list;
	while(sp!=NULL){
		if(sp->spp == n)
			return sp;
		sp = sp->next;
	}

	printf("Species %d not found\n",n);
	return NULL;

}





int stringPM::get_ecosystem(){
	l_spp *sp;
	s_ag *ag;
	int maxsp;
	int maxct;

	//reset the flags
	sp=spl->species_list;
	while(sp!=NULL){
		sp->pf=sp->anc=0;
		sp = sp->next;
	}

	//Find the current ecosystem
	ag=nowhead;
	while(ag!=NULL){
		sp=spl->species_list;
		while(sp!=NULL){
			if(sp->spp == ag->spp->spp){
				sp->pf++;
				break;
			}
			sp=sp->next;
		}
		ag=ag->next;
	}

	maxsp= spl->species_list->spp;
	maxct=spl->species_list->pf;
	sp=spl->species_list;
	while(sp!=NULL){
		if(sp->pf > maxct){
			maxsp = sp->spp;
			maxct = sp->pf;
		}
		sp=sp->next;
	}

	return maxsp;
}





/*******************************************************************************
* @brief trace the ancestry of the current population and write to a dot file
*
* @param[in] fp the file path
*
* @param[in] time the current time
*
* @param[in] step unused - TODO(sjh) refactor!
*******************************************************************************/
void stringPM::SpeciesPrintAncestryDot(FILE *fp, int time){

	l_spp *sp;
	s_ag *ag;

	const int minno = 1; //There must be more than this many individuals present for the ancestry to be plotted.....

	if(fp==NULL)
		fp=fopen("ancestry.dot","w");

	//reset the flags
	sp=spl->species_list;
	while(sp!=NULL){
		sp->pf=sp->anc=0;
		sp = sp->next;
	}

	//Find the current ecosystem
	ag=nowhead;
	while(ag!=NULL){
		sp=spl->species_list;
		while(sp!=NULL){
			if(sp->spp == ag->spp->spp){
				sp->pf++;
				break;
			}
			sp=sp->next;
		}
		ag=ag->next;
	}

	//Now trace it's origins...
	int finished = 0;
	l_spp *sp2;
	while(!finished){
		finished=1;
		sp=spl->species_list;
		while(sp!=NULL){
			if(sp->pf > minno || sp->anc){

				if(sp->pp->pp != NULL){
					sp2 = SpeciesFromListByIndex(sp->pp->pp->spp);
					if(sp2==NULL){
						printf("NULL sp2\n");fflush(stdout);
					}
					//if(sp2->sptype){//This means its not a "start" type...
						//if(!sp2->pf){
					else{	
						if(!sp2->anc){
							sp2->anc=1;
							finished=0;
						}
					}
				}
					//}
				if(sp->pp->pp != NULL){
					sp2 = SpeciesFromListByIndex(sp->pp->pa->spp);
					//if(sp2->sptype){//This means its not a "start" type...
					if(!sp2->anc){
						sp2->anc=1;
						finished=0;
					}
					//}
				}

			}
			sp=sp->next;
		}
	}

	//Now we can print it!!!


	fprintf(fp,"digraph simple_hierarchy {\noverlap = scale\n");

	//Do the timeline first:	//TIME LINE GRAPH:
	fprintf(fp,"/* TIME LINE */\n{node [shape=plaintext, fontsize=16];\n");
	fprintf(fp,"0");
	sp = spl->species_list;
	while(sp!=NULL){
		if(sp->tspp){
			if(sp->pf > minno || sp->anc){
				fprintf(fp," -> %d",sp->tspp);
			}
		}
		sp=sp->next;
	}
	//Have a line for NOW:
	fprintf(fp," -> %d",time);
	fprintf(fp,";\n}\n");


	fprintf(fp,"\n\n/* ANCESTORS: */");
	sp=spl->species_list;
	while(sp!=NULL){

		//if(sp->pf<0 || sp->pf > minno){//It's an ancestor
		if(sp->anc || sp->pf > minno){//Its an ancestor, NOT a minno
			if(!sp->pf)
				fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [style=filled,color=\"grey\",label=\"%03d\"]}\n",sp->tspp,sp->spp,sp->spp);
			else
				fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [label=\"%03d\"]}\n",sp->tspp,sp->spp,sp->spp);

			fprintf(fp,"\t\t\t cl%03d [shape=box,width=0.15,height=0.15,style=filled,color=",sp->spp);
			switch(sp->sptype){
			case 'C':
				fprintf(fp,"green");
				break;
			case 'A':
				fprintf(fp,"red");
				break;
			case 'P':
				fprintf(fp,"black");
				break;
			default:
				fprintf(fp,"grey");
			}
			fprintf(fp,",label=\"\"]\n");
			if(sp->pp->pa!=NULL){
				fprintf(fp,"\t\t\t anc%03d -> cl%03d\n",sp->pp->pa->spp,sp->spp);
				sp2=SpeciesFromListByIndex(sp->pp->pa->spp);
				if(sp2->pf>0 && sp2->pf <=minno){
					fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [label=\"%03d\"]}\n",sp2->tspp,sp2->spp,sp2->spp);
				}
			}
			if(sp->pp->pp!=NULL){
				fprintf(fp,"\t\t\t anc%03d -> cl%03d [color=grey]\n",sp->pp->pp->spp,sp->spp);
				sp2=SpeciesFromListByIndex(sp->pp->pp->spp);
				if(sp2->pf>0 && sp2->pf <=minno){
					fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [label=\"%03d\"]}\n",sp2->tspp,sp2->spp,sp2->spp);
				}
			}
			if(sp->pp->pa!=NULL || sp->pp->pp!=NULL)
				fprintf(fp,"\t\t\t cl%03d -> anc%03d\n",sp->spp,sp->spp);
		}
		sp=sp->next;
	}

	fprintf(fp,"\n\n/* CURRENT ECOSYSTEM: */");
	sp=spl->species_list;
	int minnoct=0;
	while(sp!=NULL){
		if(sp->pf>minno){
			fprintf(fp,"\t\t\t {rank=same; %d;  spp%03d [style=filled,color=black,fontcolor=white,label=\"%03d\"]}\n",time,sp->spp,sp->spp);
			fprintf(fp,"\t\t\t anc%03d -> spp%03d [color=red,headlabel=\"%d\"]\n\n",sp->spp,sp->spp,sp->pf);
		}
		//if(sp->pf>0 && sp->pf <= minno){
		//	minnoct += sp->pf;
		if(sp->pf && sp->pf <=minno){
			minnoct+= sp->pf;
		}
		sp=sp->next;
	}
	if(minnoct){
		fprintf(fp,"\t\t\t {rank=same; %d;  minnos [style=filled,color=black,fontcolor=red,label=\"rare\"]}\n",time);
		fprintf(fp,"ancminno [style=invis]");
		fprintf(fp,"\t\t\t ancminno -> minnos [color=red,headlabel=\"%d\"]\n\n",minnoct);
	}
	fprintf(fp,";\n}\n");
	fflush(fp);
	fclose(fp);

}





/*******************************************************************************
* @brief print the strings and associated data of species in the species list
*
* @param[in] fp the file path
*
* @param[in] time the current time
*
* @param[in] step unused - TODO(sjh) refactor!
*******************************************************************************/
void stringPM::SpeciesPrintStrings(FILE *fp){


	if(fp==NULL)
		fp=fopen("species.txt","w");

	l_spp *sp;


	fprintf(fp,"SppNo\tTIME\tNumber\tP1\tP2\tbtype\tlen\tSTRING\n");
	for(sp=spl->species_list;sp!=NULL;sp=sp->next){
		int len = strlen(sp->S);
		//TODO: record all parent pairs.
		s_parent *ppp;
		for(ppp = sp->pp; ppp!= NULL; ppp=ppp->next){
			if(sp->pp->pa!=NULL || sp->pp->pp!=NULL){
				fprintf(fp,"%05d\t%06d\t%04d\t%d\t%d\t%d\t%d\t%s\t%p\n",
						sp->spp, sp->tspp, sp->pf, sp->pp->pp->spp,
						sp->pp->pa->spp, sp->sptype, len,  sp->S,  &(sp->S));
			}
		}
	}

	l_spp *ss;
	s_parent *pp;

	//For debugging - print the master list of species
	fprintf(fp,"\n\nMASTER SPECIES LIST");
	fprintf(fp,"SppNo\ttP1\tP2\tbtype\tlen\tSTRING\n");
	for(ss=spl->species_list ;ss!=NULL;ss=ss->next){
		int len = strlen(ss->S);
		for(pp=ss->pp; pp!=NULL; pp=pp->next){
			if(ss->pp->pa != NULL)
				fprintf(fp,"%05d\t%d\t%d\t%d\t%d\t%s\t%p\n",ss->spp, ss->pp->pp->spp, ss->pp->pa->spp, ss->sptype, len,  ss->S,  &(ss->S));
			else
				fprintf(fp,"%05d\t_\t_\t%d\t%d\t%s\t%p\n",ss->spp, ss->sptype, len,  ss->S,  &(ss->S));
		}
	}

	fflush(fp);
	fclose(fp);
}




/*
//Related to early runs where an 'epoch' was a period of dominance of one species
void stringPM::set_epochs(){

	lastepoch=get_ecosystem();
	thisepoch=lastepoch;
	nepochs=1;
}*/





int stringPM::share_agents(s_ag **hp){

	s_ag *pa,**head,*tmp,*aa,*tmphead;
	int ntot=0,herect=0,therect=0;

	tmphead = *hp;
	pa = tmphead;
	*hp = NULL;
	aa = NULL;

	while(pa !=NULL){
		ntot++;
		float rno;

		//Work out where it's going:
		if((rno = RandomBetween0And1()) < 0.5){
			head = hp;
			therect++;
		}
		else{
			head = &nowhead;
			herect++;
		}

		//Move it and any partners.
		tmp=pa->next;
		AgentExtract(&tmphead,pa);
		AgentAppend(head,pa);

		switch(pa->status){
		case B_ACTIVE:
			aa = pa->pass;
			break;
		case B_PASSIVE:
			aa = pa->exec;
			break;
		case B_UNBOUND:
			aa = NULL;
			break;
		}

		if(aa){
			ntot++;
			if(tmp==aa)
				tmp=tmp->next;
			AgentExtract(&tmphead,aa);
			AgentAppend(head,aa);
		}

		pa=tmp;
	}

	printf("Total %d. Moved %d. Left %d\n",ntot,herect,therect);
	return 1;
}





// todo(sjh): copy_agents function in version 0.2.4 if needed





void stringPM::print_agent_cfg(FILE *fp, s_ag *pa, const int pass_index = 0){

	if(pa->status == B_ACTIVE){
		fprintf(fp,"###############  %d %s\n",pa->spp->spp,pa->spp->S);
		fprintf(fp,"REACTION_ACTIVE  %d %s irwf: ",pa->spp->spp,pa->S);
		//IRWF
		fprintf(fp,"%d %d %d ",pa->it , (int) (pa->i[0]-&(pa->pass->S[0])) , (int) (pa->i[1]-&(pa->S[0])) );
		fprintf(fp,"%d %d %d ",pa->rt , (int) (pa->r[0]-&(pa->pass->S[0])) , (int) (pa->r[1]-&(pa->S[0])) );
		fprintf(fp,"%d %d %d ",pa->wt , (int) (pa->w[0]-&(pa->pass->S[0])) , (int) (pa->w[1]-&(pa->S[0])) );
		fprintf(fp,"%d %d %d" ,pa->ft , (int) (pa->f[0]-&(pa->pass->S[0])) , (int) (pa->f[1]-&(pa->S[0])) );

		if(grid){
			fprintf(fp," grid: %d %d ",pa->x,pa->y);
		}

		fprintf(fp," pass_index: %d",pass_index);

		fprintf(fp,"\n");
		fflush(fp);
		fflush(fp);//doing it twice so we can set a breakpoint after flushing
	}
	else{
		fprintf(fp,"REACTION_PASSIVE %d %s",pa->spp->spp,pa->S);
		if(grid){
			fprintf(fp, " grid: %d %d",pa->x, pa->y);
		}

		fprintf(fp, "\n");

		fprintf(fp,"###############  %d %s\n\n",pa->spp->spp,pa->spp->S);
	}

}



void stringPM::write_extant_spp(FILE *fp){
	s_ag *pag,*bag,**agarray;
	int nag,*done,
		nreactions=0;
	bool found_pass;


	//int spc,count;
	//int finished = 0;
	//int i,found;

	nag = AgentsCount(nowhead,-1);

	fprintf(fp,"%%%%%%The number of reactants in the system:\nNUMAGENTS %d\n",nag);

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));


	//Create the array from the linked list:
	agarray = (s_ag **) malloc(nag * sizeof(s_ag *));
	pag=nowhead;
	for(int aa=0;aa<nag;aa++){
		agarray[aa] = pag;
		pag = pag->next;
	}

	for(int aa=0;aa<nag;aa++){

		pag = agarray[aa];
		int bb = -1;

		switch(pag->status){
		case B_UNBOUND:
			//Write this as a traditional agent...
			fprintf(fp,"AGENT  %s 1 Q\n",pag->S);fflush(fp);
			fprintf(fp,"GRIDPOS  %d %d\n\n",pag->x,pag->y);

			done[aa] = 1;

			break;

		case B_ACTIVE:
			bag = pag->pass;
			found_pass = false;
			for(bb = 0; bb<nag; bb++){
				if(agarray[bb] == bag){
					found_pass = true;
					break;
				}
			}
			if(!found_pass){
				printf("ERROR: unable to find passive partner\n");fflush(stdout);
			}

			//Write the index of the passive partner - to be recovered when loaded...

			fprintf(fp,"%%%%%% REACTION %d\nREACTION\n",++nreactions);
			print_agent_cfg(fp, pag, bb);


			break;
		case B_PASSIVE:


			//TODO: check this works...
			fprintf(fp,"\nPASSIVE\n");
			print_agent_cfg(fp, pag);

			break;
		}
	}

	//Tidy up..
	free(done);
	free(agarray);

}


/*******************************************************************************************/
/* THIS IS THE OLD WAY OF WRITING REACTIONS, (DON'T DELETE - THIS WAS USED FOR EVOEVO PAPER
 * WITH PAULIEN HOGEWEG

	// Now go through NOWHEAD and do each thing at a time
	// Unbound first (keeps things tidy)

	pag = nowhead;
	while(pag != NULL){

		s_ag *bag;
		int count;
		switch(pag->status){
		case B_UNBOUND:
			count = 0;

			//count any unbounds with the same sppno
			for(bag=pag;bag!=NULL;bag=bag->next){
				if(bag->spp->spp == pag->spp->spp && bag->status == B_UNBOUND){
					count++;
				}
			}
			//Write this as a traditional agent...
			fprintf(fp,"AGENT  %s %d Q\n",pag->S,count);fflush(fp);

			//record the grid position if any
			if(grid){
				for(bag=pag;bag!=NULL;bag=bag->next){
					if(bag->spp->spp == pag->spp->spp && bag->status == B_UNBOUND){
						fprintf(fp,"GRIDPOS  %d %d\n",bag->x,bag->y);fflush(fp);
					}
				}
			}

			s_ag *cag;
			//move these to nexthead
			for(bag=pag;bag!=NULL;bag=cag){
				cag = bag->next; //TODO: there has to be a better way to do this..?
				if(bag->spp->spp == pag->spp->spp && bag->status == B_UNBOUND){
					extract_ag(&nowhead,bag);
					append_ag(&(nexthead),bag);
				}
			}

			pag = nowhead;
			break;
		case B_PASSIVE:
		case B_ACTIVE:
		default:
			pag = pag->next;
			break;
		}
	}


	// Now the reactions:

	fprintf(fp,"\n\n\n%%%%%% REACTIONS  %%%%%%\n\n");
	pag = nowhead;
	int nreactions = 0;
	while(pag != NULL){

		s_ag *bag;
		switch(pag->status){
		case B_ACTIVE:
			extract_ag(&nowhead,pag);

			bag = pag->pass;
			extract_ag(&nowhead,bag);

			fprintf(fp,"%%%%%% REACTION %d\nREACTION\n",++nreactions);
			print_agent_cfg(fp, pag);
			print_agent_cfg(fp, bag);

			append_ag(&(nexthead),pag);
			append_ag(&(nexthead),bag);
			pag = nowhead;
			break;
		default:
		case B_PASSIVE:
		case B_UNBOUND:
			pag = pag->next;
			break;
		}
	}

	if(nowhead!=NULL){
		printf("ERROR - there shouldn't be any agents left in nowhead\n");
		fflush(stdout);
	}

	//Move everything back to nowhead...
	update();

	***************************************/




int stringPM::print_conf(FILE *fp){


	if(fp==NULL)
		fp=stdout;



/* TODO: record NCONTAINERS value */


/*
	%%%CELL PARAMETERS
	CELLRAD		2500
	AGRAD		10
	ENERGY		0
	NSTEPS		1200000000
	USING	/n/staff/sjh/current/ewSTRING/SMconfigs/instr_set1.mis
*/
	fprintf(fp,"%%%%%%AUTOMATICALLY GENERATED Stringmol CONFIG FILE\n");
	fprintf(fp,"%%%%%%Generated at time:\nEXTIT  		%u\n\n",timestep);
	fprintf(fp,"%%%%%%CELL PARAMETERS\n");
	fprintf(fp,"CELLRAD		%d\n",(int) cellrad);
	fprintf(fp,"AGRAD       %d\n",(int) agrad);
	fprintf(fp,"ENERGY		%d\n",(int) energy);
	fprintf(fp,"NSTEPS		%d\n",(int) nsteps);//1200000000\n");
	fprintf(fp,"DECAY       %f\n",decayrate);



	/* Mutation rate is complicated. In the original paper, there were two rates:
	 * 		Indelrate: The rate of insertion and deletion
	 * 		Subrate: The rate of substitution.
	 * 	Later papers set these to be equal, read in by the MUTATE option. This means there are three possible
	 * 	settings:
	 * 		1: original rates - these rates are hard-coded, and set if there is no MUTATE parameter
	 * 		2: mutation rate - happens if a MUTATE value is set
	 * 		3: no mutation - happens with the line MUTATE 0
	 *
	 * 	There is also a parameter in the code called 'domut' which can be set to zero at anytime to turn mutation off.
	 *
	 * 	There is no facility to set separate indel and substitution rates at present
	 *
	 */


	if((indelrate - subrate)<FLT_EPSILON){
		if(indelrate<FLT_EPSILON)
			fprintf(fp,"MUTATE		0\n");
		else
			fprintf(fp,"MUTATE		%f\n",indelrate);
	}
	else{
		//No mutation rate needs to be set: the hard-wired alife values will be loaded.
	}

	fprintf(fp,"\n%%%%%% REPORTING PARAMETERS %%%%%%\n");
	fprintf(fp,"REPORTEVERY %d\n",(int) report_every);
	fprintf(fp,"IMAGEEVERY	%d\n",(int) image_every);


	/*TODO: Need to distinguish between 'USING' and 'SUBMAT' configurations.
	 * USING: files have the '.mis' extension
	 * SUBMAT: files have the '.mtx' extension
	 * Also need to know what to do if these are missing - do we create a separate file?
	 */
	if(strstr(swt_fn,".mtx")!=NULL){
		fprintf(fp,"SUBMAT	%s\n\n",swt_fn);
	}
	else{
		if(strstr(swt_fn,".mis")!=NULL){
			fprintf(fp,"USING	%s\n\n",swt_fn);
		}
		else{
			fprintf(fp,"\n%%%%%% Warning! no SWT substitution matrix specified\n");
		}
	}



	fprintf(fp,"%%%%%%RANDOM NUMBER SEED:\nRANDSEED  		%lu\n\n\n",randseed);

	char mt_file[128];
	FILE *fp2;
	sprintf(mt_file,"RNGstate_%u.txt",timestep);
	if((fp2 = fopen(mt_file,"w"))!=NULL){
		MersenneTwisterPrintStatusToFile(fp2);
		fclose(fp2);
	}
	else{
		printf("Failed to record RNG state to file %s\n",mt_file);
	}


	fprintf(fp,"%%%%%%RNG DATA FILE:\nRNGFILE %s\n\n",mt_file);


	fprintf(fp,"%%%%%%GRID PARAMETERS\n");
	fprintf(fp,"GRIDX		%d\n",(int) grid->gridx);
	fprintf(fp,"GRIDY       %d\n\n\n",(int) grid->gridy);




/* First write the extant species and total species numbers, so we can recreate the whole run easily */
	fprintf(fp,"%%%%%% SPECIES DATA %%%%%%\n");
	fprintf(fp,"%%%%%% Species seen so far %%%%%%\n");
	fprintf(fp,"TOTSPPCT    %ld\n",spl->spp_count);

	//Now count the species extant at this time.
	int *extant;
	extant = (int *) malloc(spl->spp_count*sizeof(int));
	memset(extant,0,spl->spp_count*sizeof(int));
	l_spp **lextant;
	lextant = (l_spp **) malloc(spl->spp_count*sizeof(l_spp *));

	l_spp *psp;
	int ss=0;

	for(ss=0;ss<spl->spp_count;ss++)
		lextant[ss]=NULL;

	ss=0;
	for(psp=spl->species_list;psp != NULL; psp = psp->next){
		extant[ss]=0;
		lextant[ss++]=psp;
	}

	s_ag *pag;

	int extct=0;
	for(pag = nowhead; pag != NULL; pag = pag->next){
		//TODO: we can't rely on pag->spp->spp for this - we'll have to generate an index...
		if(!extant[pag->spp->spp]){
			for(ss=0;ss<spl->spp_count;ss++){
				if(lextant[ss]==pag->spp){
					extant[pag->spp->spp]=1;
					extct++;
					fprintf(fp,"SPECIES     %d %s\n",pag->spp->spp,pag->spp->S);
				}
			}
		}
	}
	fprintf(fp,"%%%%%% extant_spp  %d\n\n\n\n",extct);


	write_extant_spp(fp);

	free(extant);
	free(lextant);

	return 0;

}





//TODO(sjh): see if this works!
/*void stringPM::print_grid(FILE *fp){

	for(int j = 0; j< grid->gridy; j++){
		for(int i=0;i<grid->gridx;i++){
			if(grid->grid[i][j] == NULL)
				fprintf(fp," ");
			else{
				switch(grid->grid[i][j]->status){
				case B_UNBOUND:
					fprintf(fp,"c");
					break;
				case B_ACTIVE:
					fprintf(fp,"@");
					break;
				case B_PASSIVE:
					fprintf(fp,"*");
					break;
				}
			}
		}
		fprintf(fp,"|\n");
	}
}*/
