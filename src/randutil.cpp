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
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "params.h"
#include "randutil.h"

#define USING_MT 			/* This means we are using the Mersenne twister      */
#define USING_SEED_DEVRAND 	/* This means we are using dev/random to get a seed. */

#define RNG_VERBOSE

#ifdef USING_MT
#include "mt19937-2.h"
#endif





/*
 * Obtain a seed from /dev/random - better than using clock, especially for array jobs
 * todo(sjh): this will only work if /dev/random is set up!
 */
/*******************************************************************************
* @brief initialise the random number generator from the OS
*
* @details accesses /dev/random to get a seed for the Mersenne Twister
*
* @return a random integer
*******************************************************************************/
unsigned int RandomNumberFromSystem(){

	unsigned int randomData = open("/dev/urandom", O_RDONLY);
	unsigned int sysRandomInteger;
	if(!read(randomData, &sysRandomInteger, sizeof sysRandomInteger)){
		printf("WARNING!: 0 bytes read from /dev/random");
		printf("in devrandomseed(), randutil.c\n");
	}
	// you now have a random integer!
	close(randomData);

#ifdef RNG_VERBOSE
	printf("in RandomNumberFromSystem(), seed is %%d: %d  (%%u: %u)\n",
		(int) sysRandomInteger,(unsigned int) sysRandomInteger);
#endif
	return sysRandomInteger;
}





/*******************************************************************************
* @brief initialise the random number generator with a seed or time
*
* @param[in] seed used to seed the rng. if <0, chosen by the program
*
* @return the value of the seed, however it was chosen
*******************************************************************************/
unsigned int RandomInit(int seed){

	unsigned int actualSeed;

	if(seed<0){
#ifdef USING_SEED_DEVRAND
		actualSeed = RandomNumberFromSystem();
#else
		actualSeed = time(NULL);
#endif
	}else{
		actualSeed = seed;
	}

#ifdef RNG_VERBOSE
	printf("in RandomInit(), seed is  %%d: %d  (%%u: %u)\n",(int) seed,(unsigned int) seed);
#endif

#ifdef USING_MT
	SetRNGSeed(actualSeed);
#else
	srand(ActualSeed);
#endif

	return actualSeed;
}





/*******************************************************************************
* @brief initialise the random number generator with a seed or time
*
* @details uses the Mersenne Twister algorithm (TODO: check range is [0,1)
*          TODO(sjh) check why we need this as well as the above...
*
* @param[in] seed used to seed the rng. if <0, chosen by the program
*
* @return the value of the seed, however it was chosen
*******************************************************************************/
unsigned long RandomInitLong(const unsigned long *inseed){

	unsigned long seed;
	if(inseed==NULL){
#ifdef USING_SEED_DEVRAND
		seed = RandomNumberFromSystem();
		printf("in longinitmyrand, seed is %ld", (long int) seed  );
		printf("(%lu)\n",(unsigned long int) seed);
#else
		seed = time(NULL);
#endif
	}
	else{
		seed = *inseed;
	}

#ifdef USING_MT
	SetRNGSeed(seed);
#else
	srand(seed);
#endif

	return seed;
}





/*******************************************************************************
* @brief generate a random number between 0 and 1
*
* @details uses the Mersenne Twister algorithm (TODO: check range is [0,1)
*
* @return a double between 0 and 1
*******************************************************************************/
double RandomBetween0And1(){
	double x;
#ifdef USING_MT
	x = GenerateRandomDouble();
#else
	x = (double) rand() / (double) RAND_MAX;
#endif
	return x;
}





unsigned long randint(){

	unsigned long x;

#ifdef USING_MT
	x = genrandint();
#else
	x = rand();
#endif
	return x;
}





/* Create an array of random integers between the range [min,max) */
int * randintarray(const int size,const int Min,const int max){
	int i, * array;
	array = (int *) malloc(size*sizeof(int));
	for(i=0;i<size;i++)
		array[i] = Min + floor((double)max*RandomBetween0And1());
	return array;
}


/* Create an array of random integers between the range [min,max) */
int * randboolarray(const int size){
	int i, * array;
	array = (int *) malloc(size*sizeof(int));
	for(i=0;i<size;i++){
		float v=RandomBetween0And1();
		if(v<0.5)
			array[i] = 0;
		else
			array[i] = 1;
	}
	return array;
}


/*UTILITY FUNCTION FOR RE-SEEDING ON RESTART*/
//todo(sjh): needs a good tidy up!
/*******************************************************************************
* @brief get the current state of the RNG
*
* @details for the Mersenne Twister, ths is the index of the state vector array
*          which is needed to reset the state on restart (the full state vector
*          array is also needed). Other RNG engines will need other information!
*
* @return the index on the state array (for MT) -
*******************************************************************************/
int RandomNumberGeneratorGetState(){
#ifdef USING_MT
	return MersenneTwisterGetState();
#else
	printf("NOT USING MERSENNE TWISTER - CAN'T GET MTI!!\n");
	return 0;
#endif
}

/*
void set_mti(int val){
#ifdef USING_MT
	mt_set_mti(val);
#else
	printf("NOT USING MERSENNE TWISTER - CAN'T SET MTI!!\n");
#endif
}*/


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
unsigned long RandomSeedInitFromFile(char *fn, int printrandseed){

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
