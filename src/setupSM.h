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

#ifndef SETUPSM_H_
#define SETUPSM_H_

	/*Utilities */
	//For restarted runs, we need to make sure we don't overwrite...
	void FilenameGetUnused(char *fn);


	/*Flags and parameters that exist outside the stringPM object*/
	struct runparams{
		int gaqnn;		//0 or 1: whether to use the QNN measure in comass_GA
		int randseed;	//if >= 0, random number seed, if -1, use /dev/random as the seed
		int indefinite;	//if 1, run forever
		unsigned int maxnsteps;	//if 0, run forever(?) if positive, max number of steps to run
	};

    void clearfiles( char *argv[]);
	void setupSMol(struct runparams &R, int argc, char *argv[]);
	void record_spp(stringPM *A);
	void SpeciesPrintCounts(stringPM *A, int t);

	void setmaxcode(stringPM *A, int *maxcode);
	int run_one_comass_trial(const int rr, stringPM *A, int * params, struct runparams *R);
	int run_one_AlifeXII_trial(stringPM *A);

	void setmutnet(const int * mutnet, swt *blosum);

	//count species in a containers nowhead
	float ctspp(stringPM *A, const int spp);

	//get stats for evolution cf seed community
	//float * evostats(char * Afn, stringPM *B,s_sw **spp_matches, float *class_score,float *self);
	void evostats(char * Afn, stringPM *B,s_sw **spp_matches, float *self, float *gvm);

	//Different ways of loading, dependant upon the no. of arguments
	int ParametersLoadFromMainArgs(stringPM *A, int argc, char *argv[], int verbose=0);


	/* Standard initialisation of random number seed */
	void init_randseed_config(int argc, char *argv[]);

	/* Print parameters of the trial */
	void print_params(stringPM *A, int ntrials, int nsteps);

	/* set up the popdy file for writing*/
	void PopdyInitFile(stringPM *A, bool overwrite = false);

/************************************************************/



	/* Spatial Stringmol functions */
	int GridSelectRandomMooreNeighbour(const int X, const int Y, const int Xlim, const int Ylim, int *xout, int *yout);
	int StringmolSpatial(int argc, char *argv[]);
	int StringmolSpatialConfigureFromFile(const char *fn, stringPM *A, smsprun **run, int runno);
	int TimestepIncrementSpatial(stringPM *A, smsprun *run);

	/* diagnostics for spatial stringmol */
	int StringmolSpatialAncestry(int argc, char *argv[]);
	int StringmolSpatialCommunity(int argc, char *argv[]);
	int StringmolSpatialPicsFromLogs(int argc, char *argv[]);

	void PNGEncodeAndSave(const char* filename, const std::vector<unsigned char>& image, unsigned width, unsigned height);

	enum smpic{
		smpic_spp,
		smpic_len
	};

	int GridSavePNG(stringPM *A, smpic pt);
#endif /* SETUPSM_H_ */
