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


#ifndef AGENT_H_
#define AGENT_H_

//TODO(sjh): this is a tangled mess.. can we tidy?
enum s_bind{B_UNBOUND=1,B_PASSIVE=2,B_ACTIVE=3};

//Forward declare the parent structure
struct s_parent;
struct l_spp;

struct s_ag{//THIS DEFINES AN INDIVIDUAL AGENT IN A STRINGMOL SYSTEM
	char 	*S;		//The executing string;
	//char 	*comp;	//the complement of the string
	int		label;
	int 	idx;	//the index of the agent.
	int		len;	//the length of the agent string.
	int 	ect;	//count of the number of instructions executed on bind
	int		nbind;	//number of binds
	s_bind 	status;
	char *i[2];		//instruction pointer
	char *r[2]; 	//read pointer
	char *w[2]; 	//write pointer
	char *f[2];		//flow (loop) pointer
	//Toggles: passive == 0; active == 1.
	int	it;			//intruction pointer toggle
	int rt;			//read pointer toggle
	int wt;			//write pointer toggle
	int ft;			//flow (loop) pointer toggle
	s_ag		*exec;
	s_ag		*pass;
	s_ag		*next;
	s_ag		*prev;

	//Lineages:
	l_spp 		*spp;
	s_parent 	*pp;

	//biomass
	int			biomass; //the biomass created during a reaction

	//spatial
	bool 		set;
	int 		x;
	int			y;
};


int AgentRewindDanglingPtrs(s_ag* act);
int AgentCheckZeroLengthString(s_ag* act);

s_ag * AgentMake(int label, const unsigned long int agct, const unsigned int maxl0);


s_ag * AgentMakeWithSequence(char * seq, const unsigned int label,
		const unsigned int agct, const unsigned int maxl0);


int AgentUnbind(s_ag * pag);

int AgentAppend(s_ag **list, s_ag *ag);

int AgentFree(s_ag *pag);

int SpeciesListUpdate(s_ag *p, char sptype, int add, l_spp *paspp,
		l_spp * ppspp, int mass);


#endif /* AGENT_H_ */
