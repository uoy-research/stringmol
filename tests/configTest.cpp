#include "catch.hpp"

#include "../src/error_codes.h"
#include "../src/default_config.h"

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

#define CONFIG_TEST_VERBOSE

// GLOBAL VARIABLES
static const char configFileName[] =  "../tests/config/alXII_config_for_tests.conf";

TEST_CASE("loading non-default parameters works"){

	//stringPM * test_config_settings( int argc, char *argv[], int return_SM){
	//stringPM * test_config_settings( int argc, char *argv[], int return_SM){
	/** The idea here is to report the default values of the parameters, then parse the config and report them again. */

	stringPM *A;

	A = new stringPM(NULL);
	unsigned int ntrials = 0;
	unsigned int nsteps = 0;
	
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
	ParameterReadOrDefineUnsignedInt(configFileName, "NSTEPS", &nsteps, -1, 1);
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

	A->BucketReset();
	delete A;
	
}






void compare_config(stringPM *A, stringPM *B){

	/*TODO: compare non-stringPM variables
	printf("Non-stringPM variables:\n");
	if(ntrials<0)
		printf("NTRIALS     not set - the default value would be used if needed\n");
	else
		printf("NTRIALS     %d\n",ntrials);
	*/

	//load params:
	//if(A->cellrad-B->cellrad>FLT_MIN){
	//	printf("ERROR - cellrad not saved properly\n");
	//}
	REQUIRE_THAT(A->cellrad, Catch::Matchers::WithinAbs(B->cellrad,0.05));
	
	/*
	if(A->vcellrad-B->vcellrad>FLT_MIN){
		printf("ERROR - vcellrad not saved properly\n");
	}
	if(A->energy!=B->energy){
		printf("ERROR - energy not saved properly\n");
	}
	if(A->nsteps!=B->nsteps){
		printf("ERROR - nsteps not saved properly\n");
	}
	


	if(A->indelrate-B->indelrate>FLT_MIN){
		printf("ERROR - indelrate not saved properly\n");
	}
	if(A->subrate-B->subrate>FLT_MIN){
		printf("ERROR - subrate not saved properly\n");
	}
	if(A->decayrate-B->decayrate>FLT_MIN){
		printf("ERROR - decayrate not saved properly\n");
	}
	if(A->maxl0!=B->maxl0){
		printf("ERROR - maxl0 not saved properly\n");
	}
	if(A->estep!=B->estep){
		printf("ERROR - estep not saved properly\n");
	}

	return 0;
	*/

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

	SMspp        SP_A;
	SMspp        SP_B;

	stringPM     A(&SP_A);
	stringPM     B(&SP_B);
    
	FILE *fp;
	char fn[] = "test_output.cfg";

	//Load the simulation and test that the config settings are correct
	//A = test_config_settings(argc,argv,1);
	//A = new stringPM(NULL);
	A.ParametersLoad(configFileName,0,1);
	A.AgentsLoad(configFileName,NULL,0,1);

	//Write the resulting config to file
	fp = fopen(fn,"w");
	A.print_conf(fp);
	fclose(fp);


	//Load the simulation and test that the config settings are correct
	B.ParametersLoad(fn,0,1);
	B.AgentsLoad(fn,NULL,0,1);
	//compare_config(A,B);
	//printf("csc for (A,B) is %d\n",csc);
	//REQUIRE_THAT(A->cellrad, WithinAbs(B->cellrad,0.05));
	CHECK(((A.cellrad - B.cellrad)*(A.cellrad - B.cellrad)) < 0.01);

	//Run the Trial forward
	//AgentsPrint(stdout,A.nowhead,0,A.maxl);





	/*
	//Write the resulting config to file
	fp = fopen(fn1000,"w");
	A->print_conf(fp);
	fclose(fp);


	//C = test_config_settings(argc,argv2,1);
	C = new stringPM(NULL);
	C->ParametersLoad(fn1000,0,1);
	C->AgentsLoad(fn1000,NULL,0,1);
	//compare_config(A,C);
	//printf("csc for (A,C) is %d\n",csc);

	//return csc;
	*/
}





TEST_CASE("Able to reload a run with RNG info from arbitrary point"){

	int argc2 = 3;
	char **argv2;
	char fn[] = "test_output.cfg";

	argv2 = (char **) malloc(argc2*sizeof(char *));
	for(int c=0;c<argc2;c++){
		argv2[c] = (char *)malloc(FN_LEN*sizeof(char));
		memset(argv2[c],0,FN_LEN*sizeof(char));
		//sprintf(argv2[c],"%s",argv[c]);
	}
	sprintf(argv2[2],"%s",fn);

	//run_one_AlifeXII_trial(&A);
	//
	SmPm_AlifeXII(argc2, argv2);

	//todo(sjh): finish this!

}

