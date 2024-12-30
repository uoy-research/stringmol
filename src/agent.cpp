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

#include "error_codes.h"

#include "memoryutil.h"
#include "randutil.h"

#include "alignment.h"
#include "stringmanip.h"

#include "agent.h"
#include "SMspp.h"




/******************************************************************************
 * @brief reposition pointers if they are beyond the end of string post cleave
 *
 * @param[in] act the agent
 *
 * @return 0 always
 *****************************************************************************/
int AgentRewindDanglingPtrs(s_ag* act){

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
 * @brief detect zero-length strings (after a Cleave for example)
 *
 * @param[in] act the agent
 *
 * @return 0 if no problem; 1 if zero-length active string;
 *         2 zero-length passive string
 *****************************************************************************/
int AgentCheckZeroLengthString(s_ag* act){

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





/*******************************************************************************
* @brief allocate memory and set up structure of a new agent
*
* @details Create an 'agent', which is a string 'molecule'
* NB: The string is not allocated here - done outside the function
*
* @param[in] alab an integer identifier (sometimes a char)
*
* @param[in] agct a counter - not sure if it counts agents or species!
*
* @return the agent
* *****************************************************************************/
s_ag * AgentMake(int alab, const unsigned long int agct){//, const unsigned int maxl0){

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
		ag->idx = agct;// used to be (agct)++ - do this *outside*; //TODO(sjh): check this!
		ag->nbind=0;
		ag->ect=0;
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
* @brief create an agent from sequence string
*
* @param[in] seq the sequence
*
* @param[in] label the code
*
* @param[in] agct the count
*
* @param[in] maxl0 the max string length
*
* @return the agent
*******************************************************************************/
s_ag * AgentMakeWithSequence(char * seq, const unsigned int label,
		const unsigned int agct,
		const unsigned int maxl0){


	if(strlen(seq) > maxl0){
		printf("Unable to allocate enough space for this agent: \n%s\nConsider specifying MAXL in your config\n",seq);
		exit(SEQ_LEN_ERROR);
	}else{
		s_ag * ag;

		ag = AgentMake(label,agct);//,maxl0);

		ag->S =(char *) malloc(maxl0*sizeof(char));
		memset(ag->S,0,maxl0*sizeof(char));
		strncpy(ag->S,seq,maxl0-1);//active_string));
		ag->len = strlen(ag->S);

		return ag;
	}

}





/*******************************************************************************
* @brief return one of a bound pair to its unbound state
*
* @param[in] pag the agent
*
* @return the mass
*******************************************************************************/
int AgentUnbind(s_ag * pag){

	int mass = 0;

	if(pag->status==B_ACTIVE){
		mass = pag->biomass; //todo(sjh): resolve usage of biomass
		pag->biomass = 0;
	}

	pag->status = B_UNBOUND;
	pag->pass = NULL;
	pag->exec = NULL;

	pag->ect=0;

	pag->f[0] = pag->i[0] = pag->r[0] = pag->w[0] = 0;
	pag->f[1] = pag->i[1] = pag->r[1] = pag->w[1] = 0;
	pag->ft   = pag->it   = pag->rt   = pag->wt = 0;

	return mass;
}





/******************************************************************************
 * @brief append an agent to a list
 *
 * @param[in] list the list of agents
 *
 * @param[in] ag the agent
 *
 * @return 0 always
 *****************************************************************************/
int AgentAppend(s_ag **list, s_ag *ag){
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
* @brief free memory used by an agent
*
* @param[in] pag the agent to free
*
* @return 0 always
*******************************************************************************/
int AgentFree(s_ag *pag){

	if(pag->S != NULL){
		//printf("destroying agent %d, code = %s\n",pag->idx,pag->S);
		free(pag->S);
	}

	free(pag);
	//TODO: we should set this to null at the moment we have to do it after each call to this function...!
	//pag = NULL;

	return 0;
}





/*******************************************************************************
* @brief align two agents
*
* @param[in] a1 first agent
*
* @param[in] a2 second agent
*
* @param[in] sw the alignment data
*
* @param[in] blosum the alignment table
*
* @param[in] swlist list of previous alignments
*
* @return the bind probability
*******************************************************************************/
float AgentsAlign(s_ag *a1, s_ag *a2, align *sw, swt *blosum, s_sw *swlist){

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

	//if(verbose_bind){
	//	printf("Alignment:\nm1: %d to %d\nm2: %d to %d\nscore = %f\nProb = %f = %E\n",sw->s1,sw->e1,sw->s2,sw->e2,sw->score,bprob,bprob);
	//}

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
void PointerPrintOffset(FILE *fp,const char *S,const char *p,int F, char c){
	int i,n=p-S;
	if(n<0){
		printf("Problem calculating pointer location\n");
		fflush(stdout);
	}
	for(i=0;i<n;i++)
		fprintf(fp," ");
	fprintf(fp,"%c\n",F?c-32:c);
}





/*******************************************************************************
 * @brief print the current reaction state
 *
 * @param[in] fp file pointer
 *
 * @param[in] act active molecule
 *
 * @param[in] pas passive molecule
 ******************************************************************************/
void ReactionPrintState(FILE *fp, s_ag *act, s_ag *pas, const int maxl){

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





/*******************************************************************************
* @brief print all the agents from a head (nowhead or nexthead usually)
*
* @param[in] fp file pointer (including stdout)
*
* @param[in] head - nowhead or nexthead
*
* @param[in] verbose verbose output
*
* @return an agent, or NULL if not available
******************************************************************************/
void AgentsPrint(FILE *fp, s_ag *head, bool verbose, const int maxl){

	s_ag *pag;
	pag = NULL;

	while(pag!=NULL){
		if(verbose){
			switch(pag->status){
				case B_UNBOUND:
					fprintf(fp,"Agent %6d,\texec=%4d\tnbind=%3d\tUNBOUND, %s\n",
							pag->idx,pag->ect,pag->nbind,pag->S);
					break;
				case B_ACTIVE:
					ReactionPrintState(fp,pag,pag->pass,maxl);
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
 * @brief decay agent with a probability
 *
 * @param[in] pag the agent
 *
 * @param[in] decayrate the prob of decay
 *
 * @param[in] dodecay whether to decay at all
 *
 * @return 1 if decay happens, 0 if not
 *****************************************************************************/
int AgentAttemptDecay(s_ag *pag, const float decayrate, const bool dodecay){

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
		AgentFree(pag);
		return 1;
	}
	else
		return 0;
}





/*******************************************************************************
* @brief Print the state of agent with index 'idx'
*
* @param[in] fp file pointer, can be stdout for print to screen etc. NB no error checking for file stats
*
* @param[in] detail flags the level of detail
*
* @param[in] idx the index of the agent
*
* @return 0 regardless of succes (todo: fix this)
*******************************************************************************/
int AgentPrintWithIndex(FILE *fp, int detail, int idx,
		s_ag * nowhead, unsigned int maxl){
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
					ReactionPrintState(stdout,pag,pag->pass,maxl);
				else{
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tACTIVE,  %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tPASSIVE, %s\n",pag->pass->idx,pag->pass->ect,pag->pass->nbind,pag->pass->S);
				}
				//ReactionPrintState(stdout,pag,pag->pass);
				return 1;
			case B_PASSIVE:
				if(detail){
					ReactionPrintState(stdout,pag->exec,pag,maxl);
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





/*******************************************************************************
* @brief count agents
*
* @param[in] head: usually 'nowhead' or 'nexthead'
*
* @param[in] state if -1: count all; else count with a particular status
*
* @return the count
*******************************************************************************/
int AgentsCount(s_ag *head, int state){
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
s_ag * AgentSelectRandomly(s_ag *head, int state){
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





/*******************************************************************************
* @brief extract an agent from a list
*
* @param[in] list the list of agents
*
* @param[in] ag the agent
*
* @return 0 always
*******************************************************************************/
int AgentExtract(s_ag **list, s_ag *ag){

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





/*******************************************************************************
* @brief determine whether an agent is present in a list (by address)
*
* @param[in] list the list
*
* @param[in] tag the agent
*
* @return true or false
*******************************************************************************/
bool AgentAddressInList(s_ag *list,const s_ag *tag){
	s_ag *pag;

	for(pag=list;pag!=NULL;pag=pag->next){
		if(tag == pag){
			return true;
		}
	}
	return false;
}
