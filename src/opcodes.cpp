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
#include "agent.h"
#include "SMspp.h"

#include "opcodes.h"


//extern const int  maxl;

//Use this to control whether h-search is stochastic or not
#define SOFT_SEARCH

/* "PRIVATE" FUNCTIONS */
void MassTableUpdate(int * mass, const int randomOpcodeIndex,
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
*******************************************************************************/
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
*         todo(sjh): maybe some error handling here would be good!
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

	len = OpcodeTemplateLength(ip, maxl);
	tp = iptr+len;

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

	SmithWatermanAlignment(tmp,sp,&A,T,0);

	//TODO: this will always match if any symbols match. There is no stochastic element..
#ifndef SOFT_SEARCH
	if(A.match)
#else
	int l = A.e1-A.s1 < A.e2-A.s2 ? A.e1-A.s1 : A.e2-A.s2;
	
	//todo(sjh): various different permutations of bprob...!
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

	//Ensure that the toggles are set - if no match found, we currenty move
	//*F to *I - might leave it on the opposite string
	//           if it was there in the 1st place...:
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





/*******************************************************************************
* @brief execute the "move" opcode (">") in a reaction
*
* @details see technical report
*
* @param[in] act the agent
*******************************************************************************/
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





/*******************************************************************************
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
*******************************************************************************/
void OpcodeIf(s_ag * act, swt *T, const int maxl){

	char *ip,*rp;
	char tmp[maxl],tmp2[maxl];
	int i,len;
	align A;

	ip = act->i[act->it];
	rp = act->r[act->rt];
	len = OpcodeTemplateLength(ip, maxl);
	ip++;

	switch(len){

	case 0:
	case 1:
		if(!*rp)
			act->i[act->it] = ip+len+1;
		else
			act->i[act->it] = ip+len;

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
			act->i[act->it] = ip+len+1;
		else
			act->i[act->it] = ip+len;
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
*******************************************************************************/
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
//todo(sjh): should probably set safe to false if the first check fails also!
//           WRITE TESTS for this!
/*******************************************************************************
* @brief check pointer positions for a Copy operation
*
* @details First, checks if read or write head is beyond max string length;
* 			returns -1 or -2 respectively if so.
* 			Second, checks if the read pointer is 'on' a string (i.e. not
* 			pointing at \0. Variable 'Safe' set to false if so.
*
* @param[in] act pointer to the active string (from which the partner string
*            can be accessed)
*
* @return 0 if successful;
*         -1 if attempt to write beyond maxl;
*         -2 if attempt to read beyond maxl;
*******************************************************************************/
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
*          write pointer, with mutation. Has the following features:
*          1: Checks position of pointers before attempting cleave -
*          see OpcodeCopyCheckSafe. Does nothing if not "safe".
*          Handles 'granular' flag by incrementing read & write pointers
*          (or not).
*
*
*
* @param[in] act pointer to the active string (from which the partner string
*            can be accessed)
*
* @return 0 always
*******************************************************************************/
int OpcodeCopy(s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass,
		SMspp * spl, const unsigned long int timestep){//, int &finished){

	int randomOpcodeIndex;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..

	if(!domut){
		indelrate = subrate =0;
	}

	//switch(OpcodeCopyCheckSafe(act,maxl,safe)){
	//case -1:
	//	return -1;
	//case -2:
	//	return -2;
	//}
	int pposcheck = OpcodeCopyCheckSafe(act,maxl,safe);

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
				OpcodeInsertInstruction(act, (int)((float) RandomBetween0And1() * blosum->N), NULL, blosum,
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

	if(pposcheck == 0){
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
	}
	else{
		s_ag *pass;
		pass = act->pass;
		spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp,
			AgentUnbind(act),timestep,maxl+1);

		spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
			AgentUnbind(pass),timestep,maxl+1);
		//finished = 1;
	}
	return 0;
}





/*******************************************************************************
* @brief Update the "spare" masses
*
* @param[in] mass the mass array
*
* @param[in] randomOpcodeIndex - the index of the randomly selected opcode
*
* @param[in] writePtrOpcodeIndex - the index of the opcode at the write ptr
*******************************************************************************/
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




//todo(sjh): I think we need a class for OpcodeCopy, which we can then derive
//           a class for OpcodeCopyComass from
/*******************************************************************************
 * @brief Insert an opcode - incorporates comass option
 *
 * @details Often this will be a random opcode, but for testing purposes, the
 * 	        *random* part is done before calling the function
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
void OpcodeInsertInstruction(const s_ag * act, int inst_idx, int *mass, swt * blosum,
		const int writePtrOpcodeIndex){


	int randomOpcodeIndex = inst_idx;//(float) RandomBetween0And1() * blosum->N;
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
int OpcodeCopy_Comass(s_ag *act, const bool domut,float indelrate,
		float subrate, const unsigned int maxl,
		swt	*blosum, const int granular_1, long &biomass, int *mass,
		SMspp * spl, const unsigned long int timestep){

	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..

	if(!domut){
		indelrate = subrate =0;
	}

	//switch(OpcodeCopyCheckSafe(act,maxl,safe)){
	//case -1:
	//	return -1;
	//case -2:
	//	return -2;
	//}
	int pposcheck = OpcodeCopyCheckSafe(act,maxl,safe);

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
				OpcodeInsertInstruction(act, (int)((float) RandomBetween0And1() * blosum->N), mass, blosum,
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
	////update lengths
	//act->len = strlen(act->S);
	//act->pass->len = strlen(act->pass->S);
	//act->i[act->it]++;

//#ifdef VERBOSE
//	if(mut)
//	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
//#endif
//	act->biomass++;
//	biomass++;
//	return 0;


	if(pposcheck == 0){
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
	}
	else{
		s_ag *pass;
		pass = act->pass;
		spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp,
			AgentUnbind(act),timestep,maxl+1);

		spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
			AgentUnbind(pass),timestep,maxl+1);
		//finished = 1;
	}
	return 0;
}





/*
//todo(sjh): see where speig_hcopy differs from other versions of hcopy!
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
*******************************************************************************/
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





/*******************************************************************************
* @brief cleave
*
* @param[in] act the active partner agent
*
* @param[in] nexthead the agent list to append things to
*
* @param[in] spl the species list
*
* @param[in] agct the agent count
*
* @param[in] timestep the current time
*
* @param[in] maxl0 the max string length including 0 terminating char
*
* @return which agents have been destroyed (if any) 0: none; 1: active only
*         2: passive only; 3: both
*******************************************************************************/
bool OpcodeCleave(s_ag *act, s_ag *nexthead, SMspp *spl,
		unsigned long int *agct,
		const unsigned int timestep, const unsigned int maxl0){

	int destroyAction = 0,cpy;
	s_ag *c,*pass,*csite;
	bool safe_append = true;

	pass = act->pass;

	//pick the mol containing the cleave site:
	csite = act->ft?act:pass;

	if(act->f[act->ft]-csite->S < csite->len){

		//1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT

		//Can't really say what the label is easily - for ECAL, it's always pass
		c = AgentMake(pass->label, *agct++);//,  maxl0);//,1);

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
		//Check the lineage
		spl->SpeciesListUpdate(c,'C',1,act->spp,pass->spp,act->biomass,timestep,maxl0);
		act->biomass=0; //reset this; we might continue to make stuff!

		//append the agent to nexthead
		AgentAppend(&nexthead,c);

		//2: HEAL THE PARENT

		memset(act->f[act->ft],0,cpy*sizeof(char));

		csite->len = strlen(csite->S);

		if((destroyAction = AgentCheckZeroLengthString(act))){

			switch(destroyAction){
			case 1://Destroy active - only append passive
				//AgentUnbindAndSpeciesListUpdate(pass,'P',1,act->spp,pass->spp,spl);

				//found = SpeciesListUpdate(pag,sptype,update,pa,pp,mass);
				//return found;

				spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
						AgentUnbind(pass),timestep,maxl0);

				AgentAppend(&nexthead,pass);
				AgentFreeAndNull(&act);
				//act = NULL;
				safe_append = false;
				break;
			case 2://Destroy passive - only append active
				//SMAgentUnbindAndSpeciesListUpdate(act,'A',1,act->spp,pass->spp);
				spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp,
						AgentUnbind(act),timestep,maxl0);
				AgentAppend(&nexthead,act);
				AgentFreeAndNull(&pass);
				//pass = NULL;
				safe_append = false;
				break;
			case 3://Destroy both
				printf("This should never happen\n");
				//SMAgentUnbindAndSpeciesListUpdate(act,'A',1,act->spp,pass->spp);
				spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp
						,AgentUnbind(act),timestep,maxl0);
				//SMAgentUnbindAndSpeciesListUpdate(pass,'P',1,act->spp,pass->spp);
				spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
						AgentUnbind(pass),timestep,maxl0);
				AgentFreeAndNull(&act);
				//act = NULL;
				AgentFreeAndNull(&pass);
				//pass = NULL;
				safe_append = false;
				break;
			default://This can't be right can it? - NB destroyAction = 0 covered here - make explicit!
				if(act->ft == act->it){
					act->i[act->it]--;
				}
				break;
			}
		}
	}

	if(!destroyAction){//todo(sjh): this would only be called if S>maxl0..
		act->i[act->it]++;
	}

	//return destroyAction;
	//if((OpcodeCleave(act,nexthead,spl,&agct,timestep,maxl0) )){
	//	safe_append=0;
	//}
	return safe_append;
}





/*******************************************************************************
* @brief Toggle pointers
*
* @param[in] act pointer to the active string (from which the partner string
*            can be accessed)
*
* @param[in] spl the species list
*
* @param[in] timestep
*
* @param[in] maxl0 the max string length
*******************************************************************************/
void OpcodeTerminate(s_ag *act, SMspp *spl, const unsigned long int timestep,
		const unsigned int maxl0){


	s_ag *pass;
	pass = act->pass;

	spl->SpeciesListUpdate(act,'A',1,act->spp,pass->spp,
			AgentUnbind(act),timestep,maxl0);

	spl->SpeciesListUpdate(pass,'P',1,act->spp,pass->spp,
			AgentUnbind(pass),timestep,maxl0);
}
