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

#include "randutil.h"

#define USING_MT 			/* This means we are using the Mersenne twister      */
#define USING_SEED_DEVRAND 	/* This means we are using dev/random to get a seed. */

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
int RandomSeedFromSystem(){

	int randomData = open("/dev/random", O_RDONLY);
	int sjhRandomInteger;
	if(!read(randomData, &sjhRandomInteger, sizeof sjhRandomInteger)){
		printf("WARNING!: 0 bytes read from /dev/random");
		printf("in devrandomseed(), randutil.c\n");
	}
	// you now have a random integer!
	close(randomData);

	printf("in devrandomseed, seed is %d (%u)\n",
		sjhRandomInteger,(unsigned int) sjhRandomInteger);
		
	return sjhRandomInteger;
}





/*******************************************************************************
* @brief initialise the random number generator with a seed or time
*
* @details uses the Mersenne Twister algorithm (TODO: check range is [0,1)
*
* @param[in] seed used to seed the rng. if <0, chosen by the program
*
* @return the value of the seed, however it was chosen
*******************************************************************************/
int RandomInit(int seed){

	if(seed<0){
#ifdef USING_SEED_DEVRAND
		seed = RandomSeedFromSystem();
#else
		seed = time(NULL);
#endif
	}

	printf("in initmyrand, seed is %d (%u)\n",seed,(unsigned int) seed);

#ifdef USING_MT
	SetRNGSeed(seed);
#else
	srand(seed);
#endif

	return seed;
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
		seed = RandomSeedFromSystem();
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
