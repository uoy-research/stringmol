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
/* along with this program. If not, see <http://www.gnu.org/licenses/>. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "agent.h"
#include "SMspp.h"


//extern const int  maxl0;

SMspp::SMspp() {
	spp_count=1;  // Why is this 1??
	species_list = NULL;
}

SMspp::~SMspp() {
	/* nothing needs doing here at present */
}


SMspp& SMspp::operator=(const SMspp &SMspp_in){

	//SMspp SMspp_out;
	l_spp *pls,*pls2,*head2;

	head2 = NULL;
	//Straight copy of the variables
	spp_count = SMspp_in.spp_count;

	pls = SMspp_in.species_list;

	if(pls == NULL){
		species_list = NULL;
	}
	else{
		species_list = NULL;
		while(pls != NULL){
			pls2 = SpeciesMakeFromString(pls->S,pls->tspp,strlen(pls->S)+10,pls->spp);
			if(species_list==NULL){
				species_list = pls2;
				head2 = species_list;
			}
			else{
				head2->next = pls2;
				head2 = head2->next;
			}
			pls = pls->next;
		}
	}

	return *this;
}





/***********************************************
 * @brief Find parents in the parents list; make new entry if not found
 *
 * @param[in] c the species
 *
 * @param[in] paspp the active parent
 *
 * @param[in] ppspp the passive parent
 *
 * @return the created parent
 ***********************************************/
s_parent *SMspp::ParentsFindOrMake(l_spp * c, l_spp *paspp, l_spp  *ppspp){

	s_parent *pp;
	for(pp = c->pp;pp!=NULL;pp=pp->next){
		if(paspp == pp->pa)
			if(ppspp == pp->pp){
				pp->n++;
				return pp;
			}
	}

	//We'll only make it to here if a parent is not found
	pp = ParentsMake(paspp,ppspp);
	ParentsAppend(c,pp);
	return pp;
}





/***********************************************
 * @brief Append a parent a species' list of parents
 *
 * @param[in] c the species
 *
 * @param[in] pp the parent list
 *
 * @return the created parent
 ***********************************************/
void SMspp::ParentsAppend(l_spp *c, s_parent *pp){
	s_parent *oo;

	if(c->pp==NULL){
		c->pp = pp;
		return;
	}
	else{
		oo = c->pp;
		while(oo->next != NULL){
			oo=oo->next;
		}
		oo->next = pp;
	}
}






/*******************************************************************************
 * @brief create an s_parent object
 *
 * @details We have to know what the parents are to do this, and to have
 *          checked that the parent pair does not already exist.
 *
 * @param[in] paspp the active parent
 *
 * @param[in] ppspp the passive parent
 *
 * @return the created parent
 ******************************************************************************/
s_parent * SMspp::ParentsMake(l_spp * paspp, l_spp * ppspp){

	s_parent *pp;
	pp = static_cast<s_parent *>(malloc(sizeof(s_parent)));

	//if((pp->pa = getspp(paspp))==NULL)
	//		printf("OOPS! bad active parent\n");
	//if((pp->pp = getspp(ppspp))==NULL)
	//		printf("OOPS! bad passive parent\n");
	pp->pa = paspp;
	pp->pp = ppspp;
	pp->n=1;
	pp->next=NULL;

	return pp;
}





/*******************************************************************************
 * @brief Create a species struct from a character string.
 *
 * @param[in] a the stringmol agent
 *
 * @param[in] extit the timestamp
 *
 * @param[in] maxl0 the max string length + `\0`
 *
 * @return the created species
 ******************************************************************************/
l_spp * SMspp::SpeciesMakeFromString(char *S, int extit, const int maxl0, const int spno){

	//! pointer to the l_spp object
	l_spp *sp;

	//TODO: check for malloc fails here...
	sp = static_cast<l_spp *>(malloc(sizeof(l_spp)));
	sp->S=(char *) malloc(maxl0+1*sizeof(char));

	memset(sp->S,0,maxl0*sizeof(char));
	//TODO: check that this is ok ???
	strncpy(sp->S,S,maxl0);//strlen(sp->S));//l);

	sp->pp = NULL;
	//TODO: resolve what spp_count refers to
	if(spno>-1){
		//Todo - check that spno doesn't already exist:
		l_spp *lsp;
		for(lsp = species_list; lsp != NULL; lsp = lsp->next){
			if(lsp->spp == spno){
				printf("ERROR: species %d already exists! string is: %s\n",lsp->spp,lsp->S);
				return NULL; //todo handle the error more gracefully!
				//TODO: consider exit here as potential catastrophic consequences!
			}
		}
		sp->spp = spno;
		spp_count = spp_count>spno?spp_count:(spno+1);
	}
	else{
		sp->spp=spp_count++;
	}
	sp->sptype=0;
	//These are from make_lspp in the stringPM class:
	sp->count=0;
	sp->next=NULL;
	sp->tspp=extit;
	sp->biomass=0;
	
	sp->pf = 0;
	sp->anc = 0;

	return sp;
}





/***********************************************
 * @brief Create a species struct from a stringmol agent.
 *
 * @param[in] a the stringmol agent
 *
 * @param[in] extit the timestamp
 *
 * @param[in] maxl0 the max string length + `\0`
 *
 * @return the created species
 ***********************************************/
l_spp * SMspp::SpeciesMakeFromAgent(s_ag *a, int extit, const int maxl0){

	l_spp *sp;

	sp = SpeciesMakeFromString(a->S,extit,maxl0,-1);
	//Here we assume that parents of the agent are already set!
	//sp->pp = ParentsMake(a->pp->pa,a->pp->pp);
	//a->pp = sp->pp;

	a->pp = NULL;
	return sp;
}





/***********************************************
 * @brief add a species_list to the *start* of the species_list list
 *
 * @param[in] sp the species_list
 ***********************************************/
void SMspp::SpeciesPrependToList(l_spp *sp){

	sp->next = species_list;
	species_list = sp;

}





l_spp * SMspp::find_spp(char *S, const int maxl0){

	l_spp *p;
	//First, check if it's in the list:
	for(p=species_list;p!=NULL;p=p->next){
		if(!strcmp(p->S,S))
			return p;
	}

	return NULL;

}


/*gets a spp no. and makes a new one if needed */
l_spp * SMspp::getspp(s_ag *a, int extit,const int maxl0){//s_spp * paspp, s_spp * ppspp){

	l_spp *p;

	p = find_spp(a->S,maxl0);

	if(p==NULL){
		//If not found, make a new one, append it and return the address of the new.
		p = SpeciesMakeFromAgent(a,extit,maxl0);
		SpeciesPrependToList(p);
	}

	return p;
}

/*gets a spp no. and makes a new one if needed */
l_spp * SMspp::getspp_from_string(char *S, int extit,const int maxl0, const int spno){//s_spp * paspp, s_spp * ppspp){

	l_spp *p;

	p = find_spp(S,maxl0);

	if(p==NULL){
		//If not found, make a new one, append it and return the address of the new.
		p = SpeciesMakeFromString(S,extit,maxl0,spno);
		SpeciesPrependToList(p);
	}

	return p;
}





/*******************************************************************************
* @brief free the parent list
*
* @param[in] pp the parent list
*******************************************************************************/
void SMspp::ParentListFree(s_parent *pp){

	while(pp !=NULL){
		s_parent *tmp;
		tmp=pp;
		pp=pp->next;
		free(tmp);
	}
}





/*******************************************************************************
* @brief free the species object
*
* @param[in] sp the species
*******************************************************************************/
void SMspp::SpeciesFree(l_spp *sp){
	free(sp->S);
	ParentListFree(sp->pp);
	free(sp);
	//NEED ALSO TO FREE PARENT LIST!
}






/*******************************************************************************
* @brief clear the species list
*
* @return 0 always
*******************************************************************************/
int SMspp::SpeciesListClear(){
	l_spp *ps,*n;

	for(ps=species_list;ps!=NULL;){
		n = ps->next;
		SpeciesFree(ps);
		ps = n;
	}
	species_list = NULL;
	spp_count=1;  // Why is this 1??
	return 0;
}





/*******************************************************************************
* @brief print the list of species_list
*
* @param[in] fp a file pointer (can be stdout)
*
* @return the number of species_list present
*******************************************************************************/
int SMspp::SpeciesListPrint(FILE *fp){

	l_spp *ps;
	s_parent *pp;
	int count=0;

	if(fp==NULL)
		fp=stdout;

	for(ps=species_list;ps!=NULL;ps=ps->next){
			for(pp=ps->pp;pp!=NULL;pp=pp->next){
				count++;
				if(pp->pa && pp->pp)//Need to find a better way....
					fprintf(fp,"%d,%d,%d,%d,%d,%d,%s\n",ps->spp, pp->pa->spp, pp->pp->spp ,pp->n,ps->tspp,ps->biomass,ps->S);
				else{
					fprintf(fp,"%d,",ps->spp);
					if(pp->pa)
						fprintf(fp,"%d,",pp->pa->spp);
					else
						fprintf(fp,"-1,");
					if(pp->pp)
						fprintf(fp,"%d,",pp->pp->spp);
					else
						fprintf(fp,"-1,");
					fprintf(fp,"%d,%s\n",pp->n,ps->S);
				}
			}
	}
	return count;
}
