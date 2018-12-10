
//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <complex>
#include <algorithm>

using namespace std;

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

int hasTimingErrorGraph(TGraph *gr){
	Double_t *xVals = gr->GetX();
	Double_t *yVals = gr->GetX();
	int has_error=0;
	for(int i=1; i<gr->GetN(); i++){
		if(xVals[i]<xVals[i]){
			has_error=1;
			break;
		}
	}
	return has_error;
}

int isBadRun(int station, int year, int run_num){

	int found;

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
			station2_exclude.push_back(2901); //jan 17 2014 surface pulser runs //exclusion by proximity
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
		station2_exclude.push_back(4795); //accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
		station2_exclude.push_back(4797); //accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
		station2_exclude.push_back(4798); //accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
		station2_exclude.push_back(4799); //accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
		station2_exclude.push_back(4800); //accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
	
	vector <double> station3_exclude;


	if(station==2){
		found = (std::find(station2_exclude.begin(), station2_exclude.end(), run_num) != station2_exclude.end());
		// if(run_num>=3139 && run_num <=3187) found=true; //cal pulser sweep (http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014)
	}
	else if(station==3){
		found = (std::find(station3_exclude.begin(), station3_exclude.end(), run_num) != station3_exclude.end());
	}

	return found;
}