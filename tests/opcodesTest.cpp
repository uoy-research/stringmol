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
#include "../src/stringPM.h"

TEST_CASE("insertion operation works during copy"){

	//SMspp		SP;
	SMspp SP;
	stringPM A(&SP);
	const int mymaxl = 40;
	
	char a[] = "BLUBO=STRINGA";
	char b[] = "XOYHOBSTRINGB";
	
	
	s_ag *pag,*bag;
	
	pag = AgentMakeWithSequence(a,'A',1,mymaxl);
	bag = AgentMakeWithSequence(b,'B',1,mymaxl);
	
	align sw;
	
	
	//Bind the two strings
	A.AgentsAlign(pag,bag,&sw);
	
	//Check the pointers are in the right place
	
	//Set the RNG value to be the required value for the operation
	//Or call the inner operation function...
	
	//Carry out the opcodeCopy operation...
	
	
	//Teardown
	
	
}




TEST_CASE("deletion operation works during copy"){
}




TEST_CASE("substitution operation works during copy"){
}
