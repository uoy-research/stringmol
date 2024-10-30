#include <string.h>

#include "catch.hpp"

//string stuff
#include "../src/hsort.h"
#include "../src/memoryutil.h"
#include "../src/stringmanip.h"
#include "../src/params.h"
#include "../src/alignment.h"
#include "../src/instructions.h"

//metabolism stuff
#include "../src/rules.h"
#include "../src/agents_base.h"
#include "../src/SMspp.h"
#include "../src/stringPM.h"


TEST_CASE("copy operator works for SMspp class"){

	SMspp sp1,sp2;
	char *s1,*s2;
	const int mymaxl = 40;
	
	//TODO: address that this is 1 higher than it should be!
	REQUIRE(sp1.spp_count == 1);
	
	s1 = (char *) malloc(mymaxl*sizeof(char));
	s2 = (char *) malloc(mymaxl*sizeof(char));
	strcpy(s1,"STRINGA");
	strcpy(s2,"STRINGB");

	//l_spp * make_spp_from_string(char *S, int extit, const int maxl0, const int spno);	
	sp1.species = sp1.make_spp_from_string(s1,1,mymaxl,1);
	REQUIRE(sp1.spp_count == 2);
	REQUIRE(sp1.species != NULL);
	sp1.species->next = sp1.make_spp_from_string(s2,1,mymaxl,2);
	
	sp2 = sp1;
	
	REQUIRE(sp1.spp_count == 3);
	REQUIRE(sp1.spp_count == sp2.spp_count);

	REQUIRE(!strcmp(sp1.species->next->S,sp2.species->next->S));

}



TEST_CASE("copy operator works for stringPM class: copy pre-run.."){
	
	
	//SMspp		SP;
	SMspp SP;
	const int mymaxl = 40;
	char s1[mymaxl];
	
	strcpy(s1,"STRINGA");
	SP.species = SP.make_spp_from_string(s1,1,mymaxl,1);
	
	strcpy(s1,"STRINGB");
	SP.species->next = SP.make_spp_from_string(s1,1,mymaxl,1);
	
	
	stringPM	A(&SP);
	
	//Make the copy
	stringPM	B(A);


	
	REQUIRE(A.spl->spp_count == 2);
	REQUIRE(B.spl->spp_count == A.spl->spp_count);

	REQUIRE(!strcmp(A.spl->species->next->S,"BADSTRING"));//B.spl->species->next->S));
}
