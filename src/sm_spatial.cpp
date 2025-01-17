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

#include "default_config.h"
#include "error_codes.h"

#include "memoryutil.h"
#include "mt19937-2.h"
#include "randutil.h"

#include "alignment.h"
#include "agent.h"
#include "SMspp.h"

#include "rules.h"
#include "opcodes.h"
#include "stringPM.h"
#include "sm_spatial.h"

#include "lodepng.h"
#include <iostream>
#include "setupSM.h"




Stringmol_Spatial::Stringmol_Spatial(SMspp * pSP)
: stringPM(pSP)
{


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
int Stringmol_Spatial::Run(int argc, char *argv[]) {

	printf("Hello spatial stringmol world\n");

	//SMspp		SP{};
	//stringPM	A{&SP};

	smsprun *run{};
	run = NULL;

	randseed = RandomSeedInitFromFile(argv[2],0);
	ConfigureFromFile(argv[2],&run,1);

	int bt,ct{0};
	ct = AgentsCount(nowhead,-1);
	printf("Initialisation done, number of molecules is %d\n",ct);

	//This used to be called here - but better to do it before smspatial_init()
	//that way you don't reinitialize the seed after (perhaps) positioning the
	//agents on the grid
	//
	//A.randseed = init_randseed(argv[2]);


	//reload file: should have identical settings to the input conig:
	FILE *fpp{};
	char fn[128]{};
	sprintf(fn,"reload_%05u.conf",timestep);
	fpp = fopen(fn,"w");
	print_conf(fpp);
	fclose(fpp);


	//Graphic to console - sometimes useful..
	//print_grid(stdout);


//	while(nagents(nowhead,-1)){
	while((timestep < nsteps) && (ct > 0)){

		//if(!(extit%100) || extit==1){
		//if(!(extit%100)){
		if(!(timestep%image_every)){
			bt = ct - AgentsCount(nowhead,B_UNBOUND);
			printf("Step %u done, number of molecules is %d, nbound = %d\n",timestep,ct,bt);

			GridSavePNG(smpic_spp);
			GridSavePNG(smpic_len);
		}

		//When debugging you can set specific timesteps here (not elegant, but quick!)
		//if(!(A.extit%1000) || A.extit==1
		//		|| (A.extit%1000) == 999 || (A.extit%1000) == 1  || !(A.extit%97) || !(A.extit%47)
		//		|| (A.extit%1000) == 99 || (A.extit%1000) == 100|| (A.extit%1000) == 101
		//		|| (A.extit>90001 && A.extit <9000))
		if(!(timestep%report_every)){


			sprintf(fn,"out1_%05u.conf",timestep);
			fpp = fopen(fn,"w");

			print_conf(fpp);
			fclose(fpp);

			FILE *fp{};
			sprintf(fn,"splist%u.dat",timestep);
			fp = fopen(fn,"w");
			spl->SpeciesListPrint(fp);
			fclose(fp);

			spl->SpeciesPrintCounts(nowhead,popdyfn,timestep);
		}

		TimestepIncrementSpatial();

#ifdef DODEBUG
		printf("Nowhead is %p, Nexthead is %p\n",nowhead,nexthead);
		s_ag *p;
		p=nowhead;
		int mno=0;
		while(p!=NULL){
			int x,y;
			find_ag_gridpos(p,run,&x,&y);
			printf("%d, %p, [%d,%d]  status: %d, bound to %p / %p, prev = %p, next = %p\n",
			++mno,p,x,y,p->status,p->exec,p->pass,p->prev,p->next);
			p = p->next;
		}
#endif

		timestep++;
		ct = AgentsCount(nowhead,-1);

	}

	printf("FINISHED smspatial\n");
	fflush(stdout);
	return 0;
}




//todo(sjh): integrate the smsprun struct with stringPM
/*******************************************************************************
* @brief reads a config file and sets up a Spatial Stringmol run
*
* @param[in] fn the config file name
*
* @param[in] run grid data in an smsprun struct
*
* @param[in] runno the run number
*
* @return 0 always
*******************************************************************************/
int Stringmol_Spatial::ConfigureFromFile(const char *fn, smsprun **run, int runno){

	ParametersLoad(fn,0,1);
	AgentsLoad(fn,NULL,0,0);

	run_number = runno;

	//Initialize the popdy file...
	PopdyInitFile(false);

	//TODO(sjh): It's a little perverse getting this run object out, but we have to decide whether the grid is 'core' stringmol...
	*run = grid;
	//if(run == NULL){
	if(*run == NULL){
		printf("No grid data entered for spatial stringmol!\nexiting...");
		exit(SPATIAL_NO_GRID_DEFINED);
	}


	//TODO: Now we have to place each agent on the grid - use the makenext() model -
	//But we only need to do this if extit == 0, because otherwise we'll have the molecular positions...
	if(!timestep){
		while(nowhead!=NULL){
			s_ag *pag;
			pag = AgentSelectRandomly(nowhead,-1);
			AgentExtract(&(nowhead),pag);
			int found = 0;

			//Check to see if a position has been set for each molecule..
			if( pag->x > -1 ){
				if( pag->y > -1){
					AgentPlaceOnGrid(pag,*run,pag->x,pag->y);
					AgentAppend(&(nexthead),pag);
					found = 1;
				}
				else{
					printf("Bad xy position for this agent (%d,%d) \n",pag->x,pag->y);
					exit(SPATIAL_BAD_XY);
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
							AgentAppend(&(nexthead),pag);

							s_ag *bag;
							bag = AgentSelectRandomly(nowhead,-1);
							if(bag != NULL){
								AgentExtract(&(nowhead),bag);
								AgentPlaceOnGrid(bag,*run,xx,yy);
								AgentAppend(&(nexthead),bag);
							}
						}
					}
				}
			}
		}

		UpdateNowNext();
	}
	else{
		//Fill the grid
		s_ag *pag{};

		for(int x=0;x<(*run)->gridx;x++)
			for(int y=0;y<(*run)->gridy;y++)
				(*run)->grid[x][y] =NULL;

		for(pag=nowhead; pag != NULL; pag=pag->next){
			AgentPlaceOnGrid(pag,*run,pag->x,pag->y);
		}

		TimestepGridIncrement(grid);
	}


	//#ifdef DEBUG
	for(int x=0;x<(*run)->gridx;x++){
		for(int y=0;y<(*run)->gridy;y++){
			if((*run)->grid[x][y]!=NULL){
				if(!(AgentAddressInList(nowhead,(*run)->grid[x][y]))){
					printf("Agent species %d not in nowhead at %d, %d\n", (*run)->grid[x][y]->spp->spp, x, y);
				}
			}
		}
	}
	//#endif

	return 0;
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
int Stringmol_Spatial::GridSavePNG(smpic pt){

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

	//Create the PNG
	std::vector<unsigned char> image;
	image.resize(grid->gridx * grid->gridy * 4);

	int x,y,val;

	//Bit masks for the 8*8*4 rgb cube
	int rmask = 0b11100000;
	int gmask = 0b00011100;
	int bmask = 0b00000011;
	int len;


	for(x=0;x<grid->gridx;++x){
		for (y=0;y<grid->gridy;++y) {

			if(grid->grid[x][y] == NULL){

				image[4 * grid->gridx * y + 4 * x + 0] = 0;//255 * !(x & y);
				image[4 * grid->gridx * y + 4 * x + 1] = 0;//x ^ y;
				image[4 * grid->gridx * y + 4 * x + 2] = 0;//x | y;
				image[4 * grid->gridx * y + 4 * x + 3] = 255;
			}
			else{
				switch(pt){
				case smpic_len:
					len = strlen(grid->grid[x][y]->spp->S);

					if(len>69){
						image[4 * grid->gridx * y + 4 * x + 0] = 255;//255 * !(x & y);
						image[4 * grid->gridx * y + 4 * x + 1] = 255;//x ^ y;
						image[4 * grid->gridx * y + 4 * x + 2] = 255;//x | y;
						image[4 * grid->gridx * y + 4 * x + 3] = 255;

					}
					else{

						image[4 * grid->gridx * y + 4 * x + 0] = colours[len][0];//255 * !(x & y);
						image[4 * grid->gridx * y + 4 * x + 1] = colours[len][1];//x ^ y;
						image[4 * grid->gridx * y + 4 * x + 2] = colours[len][2];//x | y;
						image[4 * grid->gridx * y + 4 * x + 3] = 255;

					}
					break;
				case smpic_spp:

					val = ((grid->grid[x][y]->spp->spp)*5) % 256;

					image[4 * grid->gridx * y + 4 * x + 0] = (32  * (1+((val & rmask) >> 5)))-1;
					image[4 * grid->gridx * y + 4 * x + 1] = (32  * (1+((val & gmask) >> 2)))-1;
					image[4 * grid->gridx * y + 4 * x + 2] = (64 *  (1+((val & bmask)     )))-1;
					image[4 * grid->gridx * y + 4 * x + 3] = 255;

					break;
				}
			}
		}
	}

	char filename[128];
	switch(pt){
	case smpic_len:
		sprintf(filename,"lenframe%07u.png",timestep);
		break;
	case smpic_spp:
		sprintf(filename,"sppframe%07u.png",timestep);
		break;

	}
	PNGEncodeAndSave(filename, image, grid->gridx, grid->gridy);

	return 0;
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
int Stringmol_Spatial::TimestepIncrementSpatial(){

	//again, we follow TimestepIncrement, but are a little more careful with the binding and uncoupling
	s_ag *pag;

	//A->energy += A->estep;

	//Set the energy to the max possible number of molecules...
	//printf("Before update, energy is %d\n",A->energy);
	energy = grid->gridx * grid->gridy;
	//printf("After update, energy is %d\n",A->energy);

	int ct=0;

	while(nowhead!=NULL){

 		s_ag *bag;
		pag = AgentSelectRandomly(nowhead,-1);
		AgentExtract(&(nowhead),pag);

		//For debugging RNG diffs.
		//todo(sjh): this should be a test...?
		if(timestep == 90001){
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
			AgentExtract(&(nowhead),bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			AgentExtract(&(nowhead),bag);
			break;
		}

		if(!AgentAttemptDecaySpatial(&pag)){
			int changed = 0;
			if(energy>0){
				switch(pag->status){
				case B_UNBOUND:
					//seek binding partner, set binding states.
					//changed = A->testbind(pag);

					//int x,y;
					align sw;

					//We can only bind neighbours in the spatial model
					//find_ag_gridpos(pag,run,&x,&y);
					grid->status[pag->x][pag->y]=G_NEXT;

					if((bag = ReactionSeekRandomSpatialPartner(pag->x,pag->y))!=NULL){

						//TODO: We need to make sure that bag is in nowhead first!
						AgentExtract(&(nowhead),bag);

						//Now we've found a potential partner, we can see if it binds:
						float bprob;
						bprob = AgentsAlign(pag,bag,&sw,blosum,swlist);

						float rno;
						rno = RandomBetween0And1();
						if(rno<bprob){//Binding success!
							//figure out which is the executing string:
							ReactionSetupExecution(pag,bag,&sw);
							pag->nbind++;
							bag->nbind++;

							energy--;

							AgentAppend(&(nexthead),pag);
							AgentAppend(&(nexthead),bag);
							changed=1;
						}
					}

					break;
				case B_PASSIVE:
					//find_ag_gridpos(pag->exec,run,&x,&y);

					changed = ReactionExecuteOpcode_Spatial(pag->exec,pag);//,x,y);

					break;
				case B_ACTIVE:

					//find_ag_gridpos(pag,run,&x,&y);
					changed = ReactionExecuteOpcode_Spatial(pag,pag->pass);//,x,y);
					break;
				default:
					printf("ERROR: agent with unknown state encountered!\n");
				}
			}
			if(!changed){
				AgentAppend(&(nexthead),pag);
				if(bag!=NULL)
					AgentAppend(&(nexthead),bag);

			}
		}
	}

	UpdateNowNext();
	TimestepGridIncrement(grid);

	//#ifdef DEBUG
	for(int x=0;x<grid->gridx;x++){
		for(int y=0;y<grid->gridy;y++){
			if(grid->grid[x][y]!=NULL){
				if(!(AgentAddressInList(nowhead,grid->grid[x][y]))){
					printf("Agent species %d not in nowhead at %d, %d\n", grid->grid[x][y]->spp->spp, x, y);
				}
			}
		}
	}

	return 0;
}


int Stringmol_Spatial::ToroidalNeighbour(int coord, int iterator, int dim){
	return (coord + iterator + dim)%dim;
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
s_ag * Stringmol_Spatial::ReactionSeekRandomSpatialPartner(int x, int y){

	int i,j,xx,yy;

	//first, let's count the agents
	int count = 0;
	for(i=-1;i<2;i++){
		for(j=-1;j<2;j++){
			if( !(i == 0 && j ==0) ){
				xx = ToroidalNeighbour(x,i,grid->gridx);
				yy = (y + j + grid->gridy)%grid->gridy;
				if(grid->grid[xx][yy]!=NULL){
					if(grid->grid[xx][yy]->status == B_UNBOUND){
						if(grid->status[xx][yy] == G_NOW){
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
				xx = (x + i + grid->gridx)%grid->gridx;
				yy = (y + j + grid->gridy)%grid->gridy;
				if(grid->grid[xx][yy]!=NULL){
					if(grid->grid[xx][yy]->status == B_UNBOUND){
						if(grid->status[xx][yy] == G_NOW){
							if(count==it)
								grid->status[xx][yy] = G_NEXT;
							return grid->grid[xx][yy];
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
int Stringmol_Spatial::ReactionExecuteOpcode_Spatial(s_ag *act, s_ag *pass){//, int x, int y){

	bool safe_append=true;

	switch(*(act->i[act->it])){//*iptr[it]){

	/*************
	 *   SEARCH  *
	 *************/
	case '$':
		OpcodeSearch(act,blosum,maxl);
		break;


	/*************
	 *   MOVE  *
	 *************/
	case '>':
		OpcodeMove(act);
		break;


	/************
	 *   HCOPY  *
	 ************/
	case '='://h-copy
		OpcodeCopy(act,domut,indelrate,subrate,maxl,
						blosum,granular_1,biomass,
						spl,timestep);
		break;


	/************
	 *   INC_R  *
	 ************/
	case '+':
		OpcodeIncrementRead(act,granular_1);
		break;


	/************
	 *  TOGGLE  *
	 ************/
	case '^'://p-toggle: toggle active pointer
		OpcodeToggle(act);
		break;


	/************
	 *  IFLABEL *
	 ************/
	case '?'://If-label
			//act->i[act->it]=OpcodeIf(act->i[act->it],act->r[act->rt],act->S,A->blosum,A->maxl);
		OpcodeIf(act,blosum,maxl);
		break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			//Decide where to put the cleaved molecule
			if((/*dac = */OpcodeCleaveSpatial(act))){//,x,y))){
				//todo(sjh): Need to determine what safe_append is used for (after looking at cleave)
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

			spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp,
					AgentUnbind(act),timestep,maxl0);

			spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
					AgentUnbind(pass),timestep,maxl0);


			grid->status[act->x][act->y] = G_NEXT;
			grid->status[pass->x][pass->y] = G_NEXT;

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
		AgentAppend(&(nexthead),act);
		AgentAppend(&(nexthead),pass);
	}
	energy--;

	return 1;
}








//todo(sjh): There seem to be more than just 'spatial' differences between
//this function and AgentAttemptDecay; resolve them!
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
int Stringmol_Spatial::AgentAttemptDecaySpatial(s_ag **pag){

	float prob = decayrate;//1./pow(65,2);//4./3.); //This is now done in load_decay...

	float rno = RandomBetween0And1();

	if(rno<prob){
		//unbind_ag(pag);

		s_ag *bag;
		bag = NULL;//To prevend compiler "uninitialised" warning
		switch((*pag)->status){
		case B_UNBOUND:
			bag = NULL;
			break;
		case B_ACTIVE:
			bag = (*pag)->pass;
			break;
		case B_PASSIVE:
			bag = (*pag)->exec;
			break;
		}
		//int x,y;

		//find_ag_gridpos(pag,run,&x,&y);
		grid->grid[(*pag)->x][(*pag)->y]=NULL;
		grid->status[(*pag)->x][(*pag)->y]=G_EMPTY;

		AgentFreeAndNull(pag);
		//TODO: sort this null-ing of free'd agents out!
		//pag = NULL;

		if(bag!=NULL){

			//find_ag_gridpos(bag,run,&x,&y);
			grid->grid[bag->x][bag->y]=NULL;
			grid->status[bag->x][bag->y]=G_EMPTY;

			AgentFreeAndNull(&bag);
			//bag = NULL;
		}

		return 1;
	}
	else
		return 0;
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
int Stringmol_Spatial::OpcodeCleaveSpatial(s_ag *act){//, int x, int y){

	int dac = 0,cpy;
	s_ag *c,*pass,*csite;
	c = NULL;
	pass = act->pass;

	//pick the mol containing the cleave site:
	csite = act->ft?act:pass;

	if(act->f[act->ft]-csite->S < csite->len){

		//1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT
		c = AgentMake(pass->label,(agct)++);//,A->maxl0);//,1);

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

		cpy -= act->f[act->ft]-cs;

		if(!cpy){
			printf("ERROR: Zero length molecule definitely being created!\nbail..\n");
			AgentFreeAndNull(&c);
			//c=NULL;
		}
		else{

			strncpy(c->S,act->f[act->ft],cpy);
			c->len = strlen(c->S);

#ifdef VERBOSE
		printf("String %d created:\n%s\n",c->idx,c->S);
#endif

			//Check the lineage
			spl->SpeciesListUpdate(c,'C',1,act->spp,pass->spp,act->biomass,timestep,maxl0);
			act->biomass=0; //reset this; we might continue to make stuff!

			//TODO: place the new agent on the grid
			if((AgentPlaceInMooreNeighbourhood(c,act->x,act->y))!=NULL){//,x,y))!=NULL){
				//append the agent to nexthead
				AgentAppend(&(nexthead),c);
			}
			else{
				AgentFreeAndNull(&c);
				//c=NULL;
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
		if((dac = AgentCheckZeroLengthString(act))){
			//int x,y;
			switch(dac){
			case 1://Destroy active - only append passive
				spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
						AgentUnbind(pass),timestep,maxl0);


				AgentAppend(&(nexthead),pass);
				//find_ag_gridpos(pass,run,&x,&y);
				//run->status[x][y]=G_NEXT;
				grid->status[pass->x][pass->y]=G_NEXT;

				//find_ag_gridpos(act,run,&x,&y);
				grid->grid[act->x][act->y]=NULL;
				grid->status[act->x][act->y]=G_EMPTY;

				AgentFreeAndNull(&act);
				//act = NULL;

				break;
			case 2://Destroy passive - only append active
				spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp,
						AgentUnbind(act),timestep,maxl0);

				AgentAppend(&(nexthead),act);
				//find_ag_gridpos(act,run,&x,&y);
				grid->status[act->x][act->y]=G_NEXT;


				//find_ag_gridpos(pass,run,&x,&y);
				grid->grid[pass->x][pass->y]=NULL;
				grid->status[pass->x][pass->y]=G_EMPTY;

				AgentFreeAndNull(&pass);
				//pass = NULL;

				break;
			case 3://Destroy both
				printf("Destroying both parents after cleave!\nThis should never happen!\n");
				spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp,
						AgentUnbind(act),timestep,maxl0);
				spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
						AgentUnbind(pass),timestep,maxl0);

				AgentFreeAndNull(&act);
				//act = NULL;
				AgentFreeAndNull(&pass);
				//pass = NULL;
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
s_ag * Stringmol_Spatial::AgentPlaceInMooreNeighbourhood(s_ag *c,int x,int y){


	int xx,yy;
	int nvacant=0;
	for(int i=-1;i<2;i++){
		for(int j=-1;j<2;j++){
			xx = (x+i+grid->gridx)%grid->gridx;
			yy = (y+j+grid->gridy)%grid->gridy;
			//No need to ingnore x,y because it is occupied by the parent!
			if(grid->grid[xx][yy]==NULL)
				nvacant++;
		}
	}
	if(nvacant){
		//Decide where to put it:
		int pos = nvacant * RandomBetween0And1();
		int here=0;
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				xx = (x+i+grid->gridx)%grid->gridx;
				yy = (y+j+grid->gridy)%grid->gridy;
				//No need to ingnore x,y because it is occupied by the parent!
				if(grid->grid[xx][yy]==NULL){
					if(here == pos){
						AgentPlaceOnGrid(c,grid,xx,yy);
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
		if(grid->grid[x][y]==NULL){
			//This should never happen....
			printf("Alert! empty parent cell!\n");
			noccupied++;
		}

		int pos = nvacant * rand0to1();
		int here=0;
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				xx = (x+i+grid->gridx)%grid->gridx;
				yy = (y+j+grid->gridy)%grid->gridy;
				//No need to ingnore x,y because it is occupied by the parent!
				if(i!=0 && j!=0){
					if(grid->grid[xx][yy]!=NULL){
						if(here == pos){
							//Remove the incumbent
							grid->grid[xx][yy]=c;
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
