#include <string.h>
#include <stdio.h>

#include "catch.hpp"

#include "../src/error_codes.h"
#include "../src/default_config.h"

//string stuff
//#include "../src/hsort.h"
//#include "../src/memoryutil.h"
//#include "../src/stringmanip.h"
//#include "../src/params.h"

//metabolism stuff
#include "../src/alignment.h"
#include "../src/rules.h"
#include "../src/agent.h"
#include "../src/SMspp.h"
#include "../src/opcodes.h"
#include "../src/stringPM.h"


s_ag * SetupReactionFromStrings(char s1[], char s2[], stringPM * A){


	const int mymaxl = 200;//todo(sjh) check/warn if s1 or s2 are longer...
	s_ag *pag,*bag;

	pag = AgentMakeWithSequence(s1,'A',1,mymaxl);
	bag = AgentMakeWithSequence(s2,'B',2,mymaxl);

	align sw;

	A->spl->SpeciesListUpdate(pag,'I',1,NULL,NULL,0,A->timestep,A->maxl0);
	A->spl->SpeciesListUpdate(bag,'I',1,NULL,NULL,0,A->timestep,A->maxl0);

	//TODO: load blosum info
	A->blosum = default_table();

	//Bind the two strings
	AgentsAlign(pag,bag,&sw,A->blosum,A->swlist);

	//Check the pointers are in the right place
	A->ReactionSetupExecution(pag,bag,&sw);

	return pag;

}


TEST_CASE("insertion operation works during copy"){

	//SMspp		SP;
	SMspp SP;
	stringPM A(&SP);
	
	const int InsertedInstructionIndex = 17;

	char a[] = "X=BLUBO=STRINGA";
	char b[] = "=OYHOBZZZ";
	
	s_ag *pag;
	
	pag = SetupReactionFromStrings(a,b,&A);

	//todo(sjh): this doesn't feel like testing cos I'm writing it years
	//after, but it will help as we refactor..

	//toggle write pointer
	pag->wt = 0;
	//move to the end of the string:
	while(*(pag->w[pag->wt]))
		pag->w[pag->wt]++;

	//useful to print state while we are checking functionality
	//A.ReactionPrintState(stdout, pag, pag->pass);
	//printf("----------inserting------------------\n");

	OpcodeInsertInstruction(pag, InsertedInstructionIndex, NULL, A.blosum, -1);

	//useful to print state while we are checking functionality
	//A.ReactionPrintState(stdout, pag, pag->pass);
	
	REQUIRE(pag->pass->S[9]=='N');
}




TEST_CASE("deletion operation works during copy"){
}




TEST_CASE("copy with the write head at the length limit leaves the string terminated"){

	SMspp SP;
	stringPM A(&SP);

	const unsigned int maxl = 20;	// small limit; the test buffers are 200 long
	char a[] = "X=BLUBO=STRINGAAAAAA";	// exactly maxl symbols
	char b[] = "=OYHOBZZZ";
	REQUIRE(strlen(a) == maxl);

	s_ag *pag = SetupReactionFromStrings(a,b,&A);
	s_ag *act = (pag->status == B_ACTIVE) ? pag : pag->exec;
	REQUIRE(act != NULL);
	REQUIRE(strlen(act->S) == maxl);

	// write head at the limit on the active string, read head at its start
	act->wt = 1;
	act->w[1] = act->S + maxl;
	act->rt = 1;
	act->r[1] = act->S;

	long biomass = 0;
	OpcodeCopy(act, false, 0, 0, maxl, A.blosum, 0, biomass, A.spl, 0);

	// nothing is written at S[maxl], so the string keeps its terminator and length,
	// and the reaction ends
	REQUIRE(act->S[maxl] == 0);
	REQUIRE(strlen(act->S) == maxl);
	REQUIRE(act->status == B_UNBOUND);
}




TEST_CASE("substitution operation works during copy"){
}
