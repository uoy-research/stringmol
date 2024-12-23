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

#include "hsort.h"





/******************************************************************************
* @brief swap component of heapsort
*****************************************************************************/
void SwapByIndex(int idx[], int a, int b){
	int tmp = idx[a];
	idx[a]=idx[b];
	idx[b]=tmp;
}





/******************************************************************************
* @brief sift component of heapsort
*****************************************************************************/
void HeapsortSiftByIndex(const int a[], int idx[],int s, const int N){
	int r = s;		//	var int root := start, child
	//int c;

	while(((r*2)+1)<N){	//	while root * 2 + 1 < count {
		int c=(r*2)+1;
		if(	(c<N-1)  &&
			(a[idx[c]] < a[idx[c+1]])
			){
			c++;
		}
		if(a[idx[r]] < a[idx[c]]){
			SwapByIndex(idx,r,c);
			r = c;
		}
		else
			return;
	}
}





/******************************************************************************
* @brief heapsort one array by index value in the other
*
* @param[in] a the array to be sorted
*
* @param[in] idx the index array
*
* @param[in] N the array size
*****************************************************************************/
void HeapsortByIndex(int a[], int idx[], const int N){
	int
		s = (N/2)-1,   	//	var int   start := count ÷ 2 - 1,
		e = N-1;    	//				end := count - 1

	//"HEAPIFY" in some descriptions:
	while(s>=0){    	//	while start ≥ 0
		HeapsortSiftByIndex(a,idx,s--,N);	//		sift(a, start, count)
					//		start := start - 1
	}
	while(e>0){
		SwapByIndex(idx,e,0);
	    HeapsortSiftByIndex(a, idx, 0, e--);
	}
}


/* TODO: might be a useful diagnostic in tests...
void print_hsort_data(int a[], int idx[], const int N, FILE *out){
	int i;
	printf("i\tidx[i]\ta[idx[i]]\ta[i]\n");
	for(i=0;i<N;i++){
		printf("%d\t%d\t%d\t\t%d\n",i,idx[i],a[idx[i]],a[i]);
	}
}*/


