#include <string.h>
#include <stdio.h>

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
	sp1.species_list = sp1.SpeciesMakeFromString(s1,1,mymaxl,-1);
	REQUIRE(sp1.spp_count == 2);
	REQUIRE(sp1.species_list != NULL);
	sp1.species_list->next = sp1.SpeciesMakeFromString(s2,1,mymaxl,-1);
	
	sp2 = sp1;
	
	REQUIRE(sp1.spp_count == 3);
	REQUIRE(sp1.spp_count == sp2.spp_count);

	REQUIRE(!strcmp(sp1.species_list->next->S,sp2.species_list->next->S));

	/*TODO: The following member variables need consideration for copying or refactoring:
	 * nowhead and nexthead - should be NULL?
	 * spp_count - not needed (another spp_count is part of spl)
	 * blosum - the s-w table
	 * swlist - list of previous alignments - check this!
	 * agct - number of molecules
	 * timestep
	 * randseed - rng seed
	 * load type - no idea
	 * biomass - possibly used in comass?
	 * bstart - no idea
	 * domut - whether mutation is allowed
	 * dodecay - whether decay is allowed
	 * lastepoch thisepoch and nepochs - tracking dominant species changes
	 * mass - for comass
	 * run_number - for multiple trials
	 * swt_fn - filename for stored smith-waterman table
	 * subrate - rate of substitution
	 * indelrate - rate of indel
	 * decayrate - rate of decay
	 * verbose_bind - <-
	 * splprint - print species list - flag? limit?
	 * report_every - frequency of writing report
	 * image_every - frequency of generating images
	 * signal - ??
	 * maxl - max molecule length
	 * maxl0 - max molecule length + 1 (for string delineation)
	 * estep - number of steps to wait until adding energy
	 * granular_1 - flag for granular option
	 * popdyfn - filename for pop dynamics data
	 * linecount - ??
	 * grid - grid for spatial stringmol
	 */

}



TEST_CASE("copy operator works for stringPM class: copy pre-run.."){
	
	
	//SMspp		SP;
	SMspp SP;
	const int mymaxl = 40;
	char s1[mymaxl];
	
	strcpy(s1,"STRINGA");
	SP.species_list = SP.SpeciesMakeFromString(s1,1,mymaxl,-1);
	
	strcpy(s1,"STRINGB");
	SP.species_list->next = SP.SpeciesMakeFromString(s1,1,mymaxl,-1);
	
	
	stringPM	A(&SP);
	
	//Make the copy
	stringPM	B(A);

	
	REQUIRE(A.spl->spp_count == 3);
	REQUIRE(B.spl->spp_count == A.spl->spp_count);

	REQUIRE(!strcmp(A.spl->species_list->next->S,B.spl->species_list->next->S));
	//REQUIRE(!strcmp(A.spl->species->next->S,"BADSTRING"));//B.spl->species->next->S));
}
