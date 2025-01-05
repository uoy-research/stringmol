#include "catch.hpp"

#include "../src/mt19937-2.h"
#include "../src/randutil.h"
#include "../src/params.h"


//metabolism stuff
#include "../src/alignment.h"
#include "../src/rules.h"
#include "../src/agent.h"
#include "../src/SMspp.h"
#include "../src/opcodes.h"
#include "../src/stringPM.h"

// Writing PNGs
#include "../src/lodepng.h"
#include <iostream>

#include "../src/setupSM.h"
#include "../src/default_config.h"

#define CONFIG_TEST_VERBOSE


TEST_CASE("loading non-default parameters works"){

	//stringPM * test_config_settings( int argc, char *argv[], int return_SM){
	//stringPM * test_config_settings( int argc, char *argv[], int return_SM){
	/** The idea here is to report the default values of the parameters, then parse the config and report them again. */

	stringPM *A;

	A = new stringPM(NULL);
	unsigned int ntrials = 0;
	unsigned int nsteps = 0;
	
	static const char configFileName[] =  "../tests/config/config_for_tests.conf";

#ifdef CONFIG_TEST_VERBOSE
	printf("\nBEFORE loading the config, params are:\n");
	print_params(A,ntrials,nsteps);
	/*
	NTRIALS     0
	NSTEPS      0
	CELLRAD     2500.000000 (vcellrad = 0.000000)
	AGRAD       10.000000
	ENERGY      0
	NSTEPS      0.000000
	BLOSUM      0 size table loaded
	MUTATE      indelrate = 0.000000; subrate = 0.000000
	DECAY       0.000000
	MAXLEN      2000, (maxl0 = 2001)
	ESTEP       20
	*/
#endif


	REQUIRE(A->dodecay == 1);
	REQUIRE(A->agct == 0);
	REQUIRE(A->spp_count == 1);
	REQUIRE(A->verbose_bind == 0);
	REQUIRE(A->maxl == STRINGPM_MAXL);
	REQUIRE(A->maxl0 == STRINGPM_MAXL0);
	REQUIRE(A->estep == 20);
	REQUIRE(A->granular_1 == false);
	REQUIRE(A->report_every == 10000);
	REQUIRE(A->image_every == 100);
	REQUIRE(A->randseed == 2008);
	//todo(sjh): sort out a proper set of params and values for ALXII.. 


	//int readordef_param_int(char *fn, const char *label, int *val, const int defaultvalue, const int verbose)
	ParameterReadOrDefineUnsignedInt(configFileName, "NTRIALS", &ntrials, 1, 1);
	int nns = ParameterReadOrDefineUnsignedInt(configFileName, "NSTEPS", &nsteps, -1, 1);
	REQUIRE(nsteps == 1000);

	//A->ConfigLoad(argv[2],NULL,0,1);
	A->ParametersLoad(configFileName,0,1);
	
	//These values are copied from the configFile..
	REQUIRE(A->energy == 0);
	REQUIRE(A->estep == 2500);
	REQUIRE(A->nsteps == 1000);
	
	//todo(sjh): make sure these are loaded by ParametersLoad 
	//           and NOT by LoadAgents...! 
	//REQUIRE(A->decayrate == 0.0005); 
	//REQUIRE(A->report_every == 50);
	//REQUIRE(A->image_every == 10);
	//REQUIRE(A->subrate == 0.0002);
	//REQUIRE(A->indelrate == 0.0002);
	//REQUIRE(A->randseed == 31324);
	
	
	
	
	
	A->AgentsLoad(configFileName,NULL,0,1);
	//if(!arg_load(A, argc, argv, 0))
	//	return NULL;

#ifdef CONFIG_TEST_VERBOSE
	printf("\n\nAFTER loading the config, params are:\n");
	print_params(A,ntrials,nsteps);
	
	if(nns==1)
		printf("NSTEPS was not specified. Simulations will run indefinitely");
	printf("..c'est ca!\n\n");
#endif

	A->BucketReset();
	delete A;
	
	REQUIRE(2 == 2);
}

