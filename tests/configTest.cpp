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


	CHECK(A->dodecay == 1);
	CHECK(A->agct == 0);
	CHECK(A->spp_count == 1);
	CHECK(A->verbose_bind == 0);
	CHECK(A->maxl == STRINGPM_MAXL);
	CHECK(A->maxl0 == STRINGPM_MAXL0);
	CHECK(A->estep == 20);
	CHECK(A->granular_1 == false);
	CHECK(A->report_every == 10000);
	CHECK(A->image_every == 100);
	CHECK(A->randseed == 2008);
	//todo(sjh): sort out a proper set of params and values for ALXII.. 


	//int readordef_param_int(char *fn, const char *label, int *val, const int defaultvalue, const int verbose)
	ParameterReadOrDefineUnsignedInt(configFileName, "NTRIALS", &ntrials, 1, 1);
	int nns = ParameterReadOrDefineUnsignedInt(configFileName, "NSTEPS", &nsteps, -1, 1);
	CHECK(nsteps == 1000);

	//A->ConfigLoad(argv[2],NULL,0,1);
	A->ParametersLoad(configFileName,0,1);
	
	//These values are copied from the configFile..
	CHECK(A->energy == 0);
	CHECK(A->estep == 2500);
	CHECK(A->nsteps == 1000);
	
	//todo(sjh): make sure these are loaded by ParametersLoad 
	//           and NOT by LoadAgents...! 
	//CHECK(A->decayrate == 0.0005); 
	//CHECK(A->report_every == 50);
	//CHECK(A->image_every == 10);
	//CHECK(A->subrate == 0.0002);
	//CHECK(A->indelrate == 0.0002);
	//CHECK(A->randseed == 31324);
	
	
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
	
}






//todo(sjh): move to catch.hpp test folder
/* Test loading and saving of configs..
 * STRATEGY:
 * 		1: Load a file with known settings - see if we've got the right number.
 * 		1a: Save the file and see if the settings are the same
 * 		2: Iterate 10,000 time steps.
 * 		3: Save state.
 * 		4: Load state in a new object
 * 		5: Compare states.
 */
TEST_CASE("Loading and saving of configs is consistent"){ 
//int test_loadsave(int argc, char *argv[]){

	/*TODO: test arguments */
	stringPM *A;
	stringPM *B;
	stringPM *C;
	const int fnlen =200;
	FILE *fp;
	char fn[] = "test_output.cfg";
	char fn1000[] = "test_output_1000.cfg";
	char **argv2;

	//Load the simulation and test that the config settings are correct
	A = test_config_settings(argc,argv,1);

	//Write the resulting config to file
	fp = fopen(fn,"w");
	A->print_conf(fp);
	fclose(fp);

	argv2 = (char **) malloc(argc*sizeof(char *));
	for(int c=0;c<argc;c++){
		argv2[c] = (char *)malloc(fnlen*sizeof(char));
		memset(argv2[c],0,fnlen*sizeof(char));
		sprintf(argv2[c],"%s",argv[c]);
	}
	sprintf(argv2[2],"%s",fn);

	//Load the simulation and test that the config settings are correct
	B = test_config_settings(argc,argv2,1);

	int csc = compare_config(A,B);
	printf("csc for (A,B) is %d\n",csc);

	//Run the Trial forward
	AgentsPrint(stdout,A->nowhead,0,A->maxl);

	run_one_AlifeXII_trial(A);

	//Write the resulting config to file
	fp = fopen(fn1000,"w");
	A->print_conf(fp);
	fclose(fp);

	sprintf(argv2[2],"%s",fn1000);


	C = test_config_settings(argc,argv2,1);


	csc = compare_config(A,C);
	printf("csc for (A,C) is %d\n",csc);


	return csc;
}

