#include <string.h>
#include <stdio.h>

#include "catch.hpp"

//string stuff
//#include "../src/hsort.h"
//#include "../src/memoryutil.h"
//#include "../src/stringmanip.h"
//#include "../src/params.h"
#include "../src/alignment.h"

//metabolism stuff
#include "../src/rules.h"
#include "../src/agents_base.h"
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

	A->SpeciesListUpdate(pag,'I',1,NULL,NULL,0);
	A->SpeciesListUpdate(bag,'I',1,NULL,NULL,0);

	//TODO: load blosum info
	A->blosum = default_table();

	//Bind the two strings
	A->AgentsAlign(pag,bag,&sw);

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




TEST_CASE("substitution operation works during copy"){
}
