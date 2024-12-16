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
 * @brief handle zero-length strings (after a Cleave for example)
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
s_ag * AgentMake(int alab, const unsigned int agct, const unsigned int maxl0){

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
		ag->nbind=0;ag->ect=0;
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







//pag = AgentMakeWithSequence("BLUBO=STRINGA",'A');
s_ag * AgentMakeWithSequence(char * seq, int label, const unsigned int agct,
		const unsigned int maxl0){

	s_ag * ag;

	ag = AgentMake(label,maxl0);

	ag->S =(char *) malloc(maxl0*sizeof(char));
	memset(ag->S,0,maxl0*sizeof(char));
	strncpy(ag->S,seq,maxl0-1);//active_string));
	ag->len = strlen(ag->S);

	return ag;

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
int AgentUnbind(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp){

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





/******************************************************************************
 * @brief free memory used by an agent
 *
 * @param[in] pag the agent to free
 *
 * @return 0 always
 *****************************************************************************/
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
int SpeciesListUpdate(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp, int mass){

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
