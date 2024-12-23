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

#include "memoryutil.h"





/***********************************************
 * @brief exit if there's been an error with malloc
 ***********************************************/
void MemoryError(){

	printf("Memory allocation error\n");
	fflush(stdout);
	getchar();
	exit(6);
}





/***********************************************
 * @brief malloc, but exit if the operation fails - this is better
 *        than SIGSEV faults that are hard to trace!
 *
 * @param[in] c the species
 *
 * @param[in] pp the parent list
 *
 * @return the created parent
 ***********************************************/
void * MallocOrExit(const int number, const int size){
	void *mem;
	if((mem = malloc(number*size))==NULL)
		MemoryError();
	return mem;
}





void	** arr2alloc(const int n1, const int n2, const int size){

	int j;
	void **aa;

	aa=(void **) malloc(n1*sizeof(void *));
	for(j=0;j<n1;j++)
		aa[j] = (void *) malloc(n2*size);

	return aa;
}





/* TODO: if we need arr2free - use this as the basis...
void  intarr2free(int **aa, const int n1, const int n2){

	int i;

	for(i=0;i<n1;i++){
		free(aa[i]);
	}
	free(aa);
}*/
