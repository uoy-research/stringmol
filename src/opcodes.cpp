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

#include "memoryutil.h"
#include "randutil.h"

//string stuff
#include "stringmanip.h"
#include "alignment.h"
#include "SMspp.h"

#include "opcodes.h"


//extern const int  maxl;

//Use this to control whether h-search is stochastic or not
#define SOFT_SEARCH

/* "PRIVATE" FUNCTIONS */
void MassTableUpdate(int * mass, const int randomOpcodeIndex,
		const int writePtrOpcodeIndex);
void OpcodeInsertRandomInstruction(const s_ag * act, int *mass, swt * blosum,
		const int writePtrOpcodeIndex);



/*******************************************************************************
 * @brief calculates the length of an opcode template
 *
 * @details see technical report section 9.1
 *
 * @param[in] ip the instruction pointer
 *
 * @param[in] maxl maximum string length
 *
 * @return the length of the template
 ******************************************************************************/
int OpcodeTemplateLength(char *ip, const int maxl){

	int len=0;
	ip = ip+1;
	while(*ip > 64 && *ip <91){
		len++;
		ip++;
	}
	if(len>maxl){
		printf("Label = %d, longer than maxl (= %d!!\n",len,maxl);
	}
	return len;
}





/*******************************************************************************
* @brief execute the "search" part of the search opcode ("$")
*
* @details see technical report
*
* @param[in] iptr the instruction pointer
*
* @param[in] sp the string that the active flow pointer is pointing at
*
* @param[in] T the Smith-Waterman object
*
* @param[in] itog toggle state of the instruction pointer
*
* @param[in] ftog toggle state of the flow pointer
*
* @param[in] maxl maximum string length
*
* @return 1 always - to indicate the reaction has changed
*         todo: maybe some error handling here would be good!
*******************************************************************************/
char * OpcodeSearchInner(char *iptr, char *sp, swt *T, const int *itog,
		int *ftog,const int maxl){

	char *ip,*tp,tmp[maxl];
	ip = iptr;
	int i,len=0;
	align A;

	memset(tmp,0,maxl*sizeof(char));

	/*NOTE: We are currently searching from the start of the string with the active flow pointer.
	Perhaps we should start AT the flow pointer, and "loop around" to the beginning of the string if
	there is no match in the first part.

	So the string:

	ABCDEFGHIJKLMNOPQRSTUVWXYZ
	          f

	would be serached as if it was:

	JKLMNOPQRSTUVWXYZABCDEFGHI

	..the best match along this line would be returned. Note it is possible that we could position <f> at the
	start of the line using this technique with a little modification. But probably better to implement a "decrement"
	operator */

	/*
	//first get the length of the string
	while(*ip > 64 && *ip <91){
		len++;
		ip++;
	}
	*/
	len = OpcodeTemplateLength(ip, maxl);
	tp = iptr+len;
	//ip=iptr+1;

	if(!len){
		//Ensure that the toggles are set:
		*ftog = *itog;
		return iptr;
	}

	memset(tmp,0,128*sizeof(char));
	strncpy(tmp,iptr+1,len);
	//generate the complement:
	for(i=0;i<len;i++)
		tmp[i] = OpcodeComplement(tmp[i]);

	//SmithWaterman(tmp,sp,&A,T,0);
	/*float bprob =*/ SmithWatermanAlignment(tmp,sp,&A,T,0);

	//TODO: this will always match if any symbols match. There is no stochastic element..
#ifndef SOFT_SEARCH
	if(A.match)
#else
	int l = A.e1-A.s1 < A.e2-A.s2 ? A.e1-A.s1 : A.e2-A.s2;
	
	//TODO: various different permutations of bprob...! 
	//if(l<=2)
	//	bprob=0;
	//else
	//	bprob = pow(A.score,l)/pow(l,l);


	float s = A.score<l-1.124? A.score : l-1.124;
	float bprob = s/(l-1.124);

	float rno = RandomBetween0And1();
	if(rno<bprob)//search success!
#endif
		return tp + A.e2 - (tp-sp);

	//Ensure that the toggles are set - if no match found, we currenty move *F to *I - might leave it on the opposite string if it was there in the 1st place...:
	*ftog = *itog;
	return tp;
}





/*******************************************************************************
* @brief execute the "search" opcode ("$") in a reaction
*
* @details see technical report
*
* @param[in] act the agent
*
* @param[in] blosum the Smith Waterman blosum matrix
*
* @param[in] maxl the max line length
*******************************************************************************/
void OpcodeSearch(s_ag *act, swt *blosum, const unsigned short int maxl){
	char *cs;
	char *tmp;

	if(act->ft)
		cs = act->S;
	else
		cs = act->pass->S;
	tmp = OpcodeSearchInner(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
	act->f[act->ft] = tmp;
	act->i[act->it]++;
}






void OpcodeMove(s_ag *act){
	char *tmp;

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
}





/******************************************************************************
* @brief execute the "if" opcode ("?")
*
* @details see technical report
*
* @param[in] ip the instruction pointer
*
* @param[in] rp the read pointer
*
* @param[in] sp the string that the active flow pointer is pointing at
*
* @param[in] T the Smith-Waterman object
*
* @param[in] maxl maximum string length
*
* @return new position of the instruction pointer, 0 if "error" (can't happen)
*****************************************************************************/
char * OpcodeIf(char *ip, char *rp, char *sp, swt *T, const int maxl){

	char tmp[maxl],tmp2[maxl];
	int i,len = OpcodeTemplateLength(ip, maxl);
	ip++;
	align A;

	switch(len){

	case 0:
		if(!*rp)
			return ip+len+1;
		return ip+len;

		break;


	case 1: //Possibly switch to look at different heads here....not implemented yet....
		if(!*rp)
			return ip+len+1;
		return ip+len;
		break;

	default:

		memset(tmp ,0,maxl*sizeof(char));
		memset(tmp2,0,maxl*sizeof(char));
		strncpy(tmp,ip,len);
		//generate the complement:
		for(i=0;i<len;i++)
			tmp[i] = OpcodeComplement(tmp[i]);

		strncpy(tmp2,rp,len);

		//SmithWaterman(tmp,tmp2,&A,T,0);
		SmithWatermanAlignment(tmp,tmp2,&A,T,0);

		if(OpcodeTemplateAligns(&A,len))
			return ip+len+1;
		return ip+len;
		break;
	}
}







/*******************************************************************************
 * @brief pointer position relative to start of string
 *
 * @param[in] pag the chemical agent
 *
 * @param[in] headtype the pointer type
 *
 * @return the pointer position
 ******************************************************************************/
int PointerPosition(s_ag *pag, char headtype){

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



//todo(sjh): Lot's of scope for refactor here...
int OpcodeCopyCheckSafe(s_ag *act, const unsigned int maxl, int & safe){

	int ppos;

	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	if( (ppos = PointerPosition(act,'w'))>=(int) maxl){

		printf("Write head out of bounds: %d\n",ppos);
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';
		act->i[act->it]++;
		return -1;
	}

	if( (ppos=PointerPosition(act,'r'))>=(int) maxl){

		printf("Read head out of bounds: %d\n",ppos);
		act->i[act->it]++;
		return -2;
	}

	// TODO(sjh): Are we handling the above errors ok?
	// Check that we aren't off the end of the string,
	// but within the allocated memory:
	if(*(act->r[act->rt]) == 0){
		safe = 0;
	}

	return 0;
}




// todo(sjh): These values should be recorded elsewhere!
// MUTATION RATES:
// THESE ARE HARD-CODED FOR NOW - THEY SHOULD BE DERIVED FROM THE BLOSUM SOMEHOW...
// const float indelrate = 0.0005,		subrate=0.375;//0.0749/2
// const float indelrate = 0.0000306125,	subrate=0.01;//0.02
// const float indelrate = 0.00006125,	subrate=0.05;//0.02
// const float indelrate = 0.000125,		subrate=0.1;//0.02
// const float indelrate = 0.000005,		subrate=0.0375;//0.0749/2
// const float indelrate = 0., 			subrate = 0.;
//
// todo(sjh): tidy up the mix of setting "safe" and returning odd values!
/*******************************************************************************
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
 ******************************************************************************/
int OpcodeCopy(s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass){

	int randomOpcodeIndex;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..

	if(!domut){
		indelrate = subrate =0;
	}

	switch(OpcodeCopyCheckSafe(act,maxl,safe)){
	case -1:
		return -1;
	case -2:
		return -2;
	}

	if(safe){






		rno=RandomBetween0And1(); //see if we are overwriting or not:
		if(rno<indelrate){//INDEL

			//should follow the blosum table for this....
			rno=RandomBetween0And1();
			if(rno<0.5){//insert

				//first do a straight copy..
				*(act->w[act->wt])=*(act->r[act->rt]);

				//increment w pointer for the insertion
				act->w[act->wt]++;

				//Do the insertion
				OpcodeInsertRandomInstruction(act, NULL, blosum,
						-1);

				//increment w again (unless we are doing granular!)
				if(granular_1==0){
					act->w[act->wt]++;
				}
			}
			else{//delete by moving the iptr and not writing anything..
				act->i[act->it]++;
			}

			if(granular_1==0){
				act->r[act->rt]++;
			}
		}
		else{
			if(rno<subrate+indelrate){//INCREMENTAL MUTATION
				
				randomOpcodeIndex = OpcodeAdjacent(*(act->r[act->rt]),blosum);
				*(act->w[act->wt])=randomOpcodeIndex;
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





void MassTableUpdate(int * mass, const int randomOpcodeIndex,
		const int writePtrOpcodeIndex){

	//if the write pointer was on the string
	if(!(writePtrOpcodeIndex<0)){
		//add to the free mass of that opcode
		mass[writePtrOpcodeIndex]++;
	}
	//decrement the free mass of the random opcode
	mass[randomOpcodeIndex]--;

}





/*******************************************************************************
 * @brief Insert a random opcode - incorporates comass option
 *
 * @details
 *
 * @param[in] act pointer to the active string (from which the partner string
 *            can be accessed)
 *
 * @param[in] mass the 'free' masses of each opcode (or NULL if not using
 *            comass)
 *
 * @param[in] randomOpcodeIndex the index of a random opcode
 *            in the blosum table
 *
 * @param[in] writePtrOpcodeIndex the index of the opcode at the write pointer
 *            in the blosum table
 *
 ******************************************************************************/
void OpcodeInsertRandomInstruction(const s_ag * act, int *mass, swt * blosum,
		const int writePtrOpcodeIndex = -1){


	int randomOpcodeIndex = (float) RandomBetween0And1() * blosum->N;
	bool update = false;

	if(mass != NULL){
		//if there's mass for this symbol:
		if(mass[randomOpcodeIndex]){
			//we'll write it to the string later
			update = true;
			MassTableUpdate(mass,randomOpcodeIndex,writePtrOpcodeIndex);
		}
	}else{
		update = true;
	}

	if(update){
		//insert the random instruction
		*(act->w[act->wt])=blosum->key[randomOpcodeIndex];
	}

}





/*******************************************************************************
 * @brief Copy operator "\=" under conservation of mass
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
 ******************************************************************************/
int OpcodeComassCopy(s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass, int *mass){

	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..

	if(!domut){
		indelrate = subrate =0;
	}

	switch(OpcodeCopyCheckSafe(act,maxl,safe)){
	case -1:
		return -1;
	case -2:
		return -2;
	}

	if(safe){

		int randomOpcodeIndex,writePtrOpcodeIndex=-1;
		if(*(act->w[act->wt])){
			writePtrOpcodeIndex=OpcodeIndex(*(act->w[act->wt]),blosum);
		}

		rno=RandomBetween0And1(); //see if we are overwriting or not:
		if(rno<indelrate){//INDEL

			//should follow the blosum table for this....
			rno=RandomBetween0And1();
			if(rno<0.5){//insert

				//first do a straight copy..
				*(act->w[act->wt])=*(act->r[act->rt]);

				//increment w pointer for the insertion
				act->w[act->wt]++;

				//Do the insertion
				OpcodeInsertRandomInstruction(act, mass, blosum,
						writePtrOpcodeIndex);

				//increment w again (unless we are doing granular!)
				if(granular_1==0){
					act->w[act->wt]++;
				}
			}
			else{//delete by moving the iptr and not writing anything..
				act->i[act->it]++;
			}

			if(granular_1==0){
				act->r[act->rt]++;
			}
		}
		else{
			if(rno<subrate+indelrate){//INCREMENTAL MUTATION

				cidx = OpcodeAdjacent(*(act->r[act->rt]),blosum);
				randomOpcodeIndex = OpcodeIndex(cidx,blosum);
				if(mass[randomOpcodeIndex]){
					*(act->w[act->wt])=cidx;
					act->w[act->wt]++;
					if(!(writePtrOpcodeIndex<0)){
						mass[writePtrOpcodeIndex]++;
					}
					mass[randomOpcodeIndex]--;
				}
				else{//
					//simply increment the read head without doing anything else
				}
				act->r[act->rt]++;//possible deletion here...
			}
			else{//NO MUTATION (but possible sub via comass effects)

				randomOpcodeIndex = OpcodeIndex(*(act->r[act->rt]),blosum);
				if(mass[randomOpcodeIndex]){
					*(act->w[act->wt])=*(act->r[act->rt]);
					act->w[act->wt]++;
					if(!(writePtrOpcodeIndex<0)){
						mass[writePtrOpcodeIndex]++;
					}
					mass[randomOpcodeIndex]--;
				}
				else{
					cidx = OpcodeAdjacent(*(act->r[act->rt]),blosum);
					randomOpcodeIndex = OpcodeIndex(cidx,blosum);
					if(mass[randomOpcodeIndex]){
						*(act->w[act->wt])=cidx;
						act->w[act->wt]++;
						if(!(writePtrOpcodeIndex<0)){
							mass[writePtrOpcodeIndex]++;
						}
						mass[randomOpcodeIndex]--;
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






/*******************************************************************************
 * @brief Increment Read operator "+"
 *
 * @param[in] act pointer to the active string (from which the partner string
 *            can be accessed)
 *
 * @param[in] granular_1 whether we are doing granular stringmol
 *
 * @return 0 if successful;
 *         -1 if attempt to write beyond maxl;
 *         -2 if attempt to read beyond maxl;
 ******************************************************************************/
void OpcodeIncrementRead(s_ag *act, bool granular_1){
	if(granular_1==1){
		char *tmp;

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
}





/*******************************************************************************
 * @brief Toggle pointers
 *
 * @param[in] act pointer to the active string (from which the partner string
 *            can be accessed)
 ******************************************************************************/
void OpcodeToggle(s_ag *act){
	char *tmp;
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
}
