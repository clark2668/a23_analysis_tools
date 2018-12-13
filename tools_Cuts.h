
//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <complex>
#include <algorithm>

#include "TGraph.h"

using namespace std;

/*
	input: tgraph of an event
	output: 0 (has no timing error), 1 (has timing error)

	function: checks to see if ever x_j < x_j+1
			which otherwise causes the interpolator to fail
			is called by hasTimingErrorMultiGraph
*/

int hasTimingErrorGraph(TGraph *gr){
	Double_t *xVals = gr->GetX();
	Double_t *yVals = gr->GetX();
	int has_error=0;
	for(int i=1; i<gr->GetN(); i++){
		if(xVals[i]<xVals[i-1]){
			has_error=1;
			break;
		}
	}
	return has_error;
}

/*
	input: vector of graphs for an event
	output: 0 (has no timing error), 1 (has timing error)

	function: checks to see if ever x_j < x_j+1
			which otherwise causes the interpolator to fail
*/

int hasTimingErrorMultiGraph(vector<TGraph*> grs){
	int event_has_error=0;
	for(int i=0; i<grs.size(); i++){
		int this_has_error = hasTimingErrorGraph(grs[i]);
		if(this_has_error==1){
			event_has_error==1;
			break;
		}
	}
	return event_has_error;
}

/*
	input: station, year, run number
	output: 0 (is good run), 1 (is bad run)

	function: looks through list of known "bad runs"
			reports if the run you gave it is bad
*/

int isBadRun(int station, int year, int run_num){

	int found;

	/*
	station 2 exclusion
	*/

	vector <double> station2_exclude;
		
		/*
		2014 surface pulsing
			originally flagged by 2884, 2895, 2903, 2912, 2916
			going to throw all runs jan 14-20 (<1 week of livetime, oh well...)
		*/
		station2_exclude.push_back(2884); //jan 14 2014 surface pulser runs //actual problem causer
			station2_exclude.push_back(2885); //jan 16 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2889); //jan 16 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2890); //jan 16 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2891); //jan 16 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2893); //jan 16 2014 surface pulser runs //exclusion by proximity
		station2_exclude.push_back(2895); //jan 16 2014 surface pulser runs //actual problem causer
			station2_exclude.push_back(2898); //jan 16 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2900); //jan 17 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2901); //jan 17 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2902); //jan 17 2014 surface pulser runs //exclusion by proximity
		station2_exclude.push_back(2903); //jan 18 2014 surface pulser runs //actual problem causer
			station2_exclude.push_back(2905); //jan 18 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2906); //jan 18 2014 surface pulser runs //exclusion by proximity
			station2_exclude.push_back(2907); //jan 18 2014 surface pulser runs //exclusion by proximity
		station2_exclude.push_back(2912); //jan 19 2014 surface pulser runs //actual problem causer
			station2_exclude.push_back(2915); //jan 18 2014 surface pulser runs //exclusion by proximity
		station2_exclude.push_back(2916); //jan 20 2014 surface pulser runs //actual problem causer
			station2_exclude.push_back(2918); //jan 20 2014 surface pulser runs

		/*
		2015 icecube deep pulsing
			4787 is the "planned" run
			4795,4797-4800 were accidental
		*/
		station2_exclude.push_back(4785); //accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001719/003/181001_ARA02AnalysisUpdate.pdf, slide 38)
		station2_exclude.push_back(4787); //deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
		for(int i=4795; i<=4800; i++){ station2_exclude.push_back(i); }

		/*
		2015 noise source tests
			January 2015
			http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2015
		*/
		for(int i=4820; i<=4825; i++){ station2_exclude.push_back(i); }
		for(int i=4850; i<=4854; i++){ station2_exclude.push_back(i); }
		for(int i=4749; i<=4936; i++){ station2_exclude.push_back(i); }
		for(int i=5210; i<=5277; i++){ station2_exclude.push_back(i); }

	/*
	station 3 exclusion
	*/
	
	vector <double> station3_exclude;
		/*
		2015 surface or deep pulsing
			got through cuts
			happened jan5-6
			waveforms clearly show double pulses or things consistent with surface pulsing
		*/
		station3_exclude.push_back(3811); //deep pulser run
			station3_exclude.push_back(3810); //elminated by proximity to deep pulser run
			station3_exclude.push_back(3820); //elminated by proximity to deep pulser run

		/*
		2015 noise source tests
			January 2015
			http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2015
		*/
		for(int i=3844; i<=3860; i++){ station3_exclude.push_back(i); }
		for(int i=3881; i<=3891; i++){ station3_exclude.push_back(i); }
		for(int i=3916; i<=3918; i++){ station3_exclude.push_back(i); }
		for(int i=3920; i<=3975; i++){ station3_exclude.push_back(i); }
		for(int i=4009; i<=4073; i++){ station3_exclude.push_back(i); }


	/*
	Check for run exclusion
	*/

	if(station==2){
		found = (std::find(station2_exclude.begin(), station2_exclude.end(), run_num) != station2_exclude.end());
	}
	else if(station==3){
		found = (std::find(station3_exclude.begin(), station3_exclude.end(), run_num) != station3_exclude.end());
	}

	return found;
}