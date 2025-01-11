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

#include "memoryutil.h"
#include "mt19937-2.h"
#include "randutil.h"

#include "alignment.h"
#include "agent.h"
#include "SMspp.h"

#include "stringPM.h"
#include "sm_spatial.h"

#include "lodepng.h"
#include <iostream>
#include "setupSM.h"








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

	SMspp		SP{};
	stringPM	A{&SP};

	smsprun *run{};
	run = NULL;

	A.randseed = RandomSeedInitFromFile(argv[2],0);
	ConfigureFromFile(argv[2],&A,&run,1);

	int bt,ct{0};
	ct = AgentsCount(A.nowhead,-1);
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
			bt = ct - AgentsCount(A.nowhead,B_UNBOUND);
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
		ct = AgentsCount(A.nowhead,-1);

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
* @param[in] A the stringPM object
*
* @param[in] run grid data in an smsprun struct
*
* @param[in] runno the run number
*
* @return 0 always
*******************************************************************************/
int Stringmol_Spatial::ConfigureFromFile(const char *fn, stringPM *A, smsprun **run, int runno){

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
			pag = AgentSelectRandomly(A->nowhead,-1);
			AgentExtract(&(A->nowhead),pag);
			int found = 0;

			//Check to see if a position has been set for each molecule..
			if( pag->x > -1 ){
				if( pag->y > -1){
					AgentPlaceOnGrid(pag,*run,pag->x,pag->y);
					AgentAppend(&(A->nexthead),pag);
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
							AgentAppend(&(A->nexthead),pag);

							s_ag *bag;
							bag = AgentSelectRandomly(A->nowhead,-1);
							if(bag != NULL){
								AgentExtract(&(A->nowhead),bag);
								AgentPlaceOnGrid(bag,*run,xx,yy);
								AgentAppend(&(A->nexthead),bag);
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
				if(!(AgentAddressInList(A->nowhead,(*run)->grid[x][y]))){
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
int Stringmol_Spatial::GridSavePNG(stringPM *A, smpic pt){


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
