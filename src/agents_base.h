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


struct s_ix{
	int n;
	float prob;
	int start;
	int stop;
	int	label;
	s_ix *next;
};


class agents_base{

	public:

		//for testing propensity
		int  *bct; //count of times we've tried for a particular reagent count
		int  *bpp; //count of times we've been closenough to run a bind...
		int  bmax; //max reagents we are going to bother with.
		void 	PropensityRecord(int N,int X);
		void 	PropensityPrint(FILE *fp);
		int 	proper_prop(const int n);
		int 	PropensityEquation(const int n);

		s_ix *ifxhead;

		float cellrad;  //The radius of the "cell"
		float agrad;    //The active radius of the agent
		float vcellrad;	//The radius as a function of the number of div count things
		float move;		//The amount an agent can move. Set to 1.1 * agrad...
		long energy;	//The initial energy present in the system
		float nsteps;	//The number of steps (SHOULD BE AN INT!)

		//diagnostics
		int *tr; 		//count of the number times a rule has been tried
		int *fr;		//count of the number times a rule has fired.

		//creators and destructors
		agents_base();
		//TODO(sjh): the destructor is virtual because I got a warning otherwise...need to check this out..?
		virtual ~agents_base();
		void ParametersSetDefaults();
		void clearout(int verbose);

		//fromfile stuff
		void ConfigLoad(const char *fn, char *fninput, int test, int verbose);
		//void load(const char *fn, int test, int verbose);

		int ParametersLoad(const char *fn, int test, int verbose);
		int load_division(char *fn);
		int load_replenish(char *fn);

		//Diagnoistics
		void print_agents_header(FILE *fp);
		void print_agents_count(FILE *fp);
		//void update_aac();

		//Counting and indexing
		int aac_count(int lab);
		int nttindex(const int label);
		void make_btab(rules *rset);
		void make_com();
		void make_dcom();
		void makefr(int nr);
		void printfr(FILE *fp, rules *rset);

		//Iteration stuff
		//virtual void TimestepIncrement(rules *rset)=0;
		virtual void TimestepIncrement()=0;

		//Cell division
		int divide_conditions(int time);

};

#endif /*AGENT_H_*/
