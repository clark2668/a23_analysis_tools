
//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <complex>
#include <algorithm>

#include "TGraph.h"
#include "TMath.h"

#include "AraGeomTool.h"

using namespace std;

/*
	input: station, config, pulser, pol, startX, stopX, startY, stopY (the last four are the bounds of the projection region for projecting the TH2D into TH1D)
	output: changes the value of thetas and phis

	function: Returns the bins for plotting the cal pulser cut
				if "inBins" is true, the answer will be reported as a bin number from 0->360/0->180
			if "inBins" is false, the answer will be reported as a number in degrees from -180->180/-90->90
*/

void getCalCutPlotBoundary(int station, int config, int pulser, int pol, bool inBins, int &startX, int &stopX, int &startY, int &stopY){
	if(station==2){
		if(pulser==5){
			if(pol==0){
				if(config==2){
					startX=-35;
					stopX=-15;
					startY=-35;
					stopY=-15;
				}
			}	
		}
		if(pulser==6){
			if(pol==0){
				if(config==1 || config==2 || config==3 || config==4){
					startX=58;
					stopX=72;
					startY=-4;
					stopY=14;
				}
			}
		}

	}
	if(inBins){
		startX+=180;
		stopX+=180;
		startY+=90;
		stopY+=90;
	}
}

/*
	input: station, config, pulser, pol, startPhi, stopPhi, startTheta, stopTheta (the last four are the bounds of the fit region for computing the cut)
	output: changes the value of thetas and phis

	function:  looks up the boundaries for the gaussian fits used to set the cal pulser cut boxes
*/


void getCalFitRange(int station, int config, int pulser, int pol, double &startPhi, double &stopPhi, double &startTheta, double &stopTheta){
	if(station==2){
		if(pulser==5){
			if(pol==0){
				if(config==2){
					startPhi=-30.;
					stopPhi=-18.;
					startTheta=-30.;
					stopTheta=-16.;
				}
			}
		}
		if(pulser==6){
			if(pol==0){
				if(config==1 || config==2 || config==3 || config==4){
					startPhi=62.;
					stopPhi=70.;
					startTheta=2.;
					stopTheta=10.;
				}
			}
		}
	}
}

/*
	input: station, pulser (5 or 6), polarization, theta and phi
	output: changes the value of theta and phi to where the cal pulser is

	function: computes the "truth" theta and phi of the cal pulser
*/

void getRealLocation(int station, int pulser, int pol, double &theta, double &phi){
	
	int cal_ant_index;
	if(pulser==5){
		if(pol==0) cal_ant_index=1;
		else if(pol==1) cal_ant_index=0;
	}
	else if(pulser==6){
		if(pol==0) cal_ant_index=3;
		else if(pol==1) cal_ant_index=4;
	}

	AraGeomTool *araGeom = AraGeomTool::Instance();

	//compute the average depth of the station
	double antenna_average[3]={0.};
	for(int i=0; i<16; i++){
		for(int ii=0; ii<3; ii++){
			antenna_average[ii]+=(araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[ii]);
		}
	}
	for(int ii=0; ii<3; ii++){
		antenna_average[ii]/=16.;
	}
	double X = araGeom->getStationInfo(station)->getCalAntennaInfo(cal_ant_index)->antLocation[0];
	double Y = araGeom->getStationInfo(station)->getCalAntennaInfo(cal_ant_index)->antLocation[1];
	double Z = araGeom->getStationInfo(station)->getCalAntennaInfo(cal_ant_index)->antLocation[2];
	phi = TMath::ATan2(Y-antenna_average[1],X-antenna_average[0])*TMath::RadToDeg();
	double depth_diff = Z-antenna_average[2];
	double horz_dist = sqrt(pow((X-antenna_average[0]),2.0)+pow((Y-antenna_average[1]),2.0));
	theta = TMath::ATan2(depth_diff,horz_dist)*TMath::RadToDeg();
}

/*
	input: vector of graphs for an event
	output: 0 (has no error), 1 (has less than 550 points error)

	function: checks to see if any waveform is < 550 samples
*/

int hasShortWaveformMultiGraph(vector <TGraph*> grs){
	int event_has_error=0;
	for(int i=0; i<grs.size(); i++){
		if(grs[i]->GetN()<550.) event_has_error++;
	}
	return event_has_error;
}

/*
	input: station, unixTime (UTC)
	output: 0 (is good time), 1 (is bad time)

	function: checks if the livetime is bad
*/

int isBadLivetime(int station, int unixTime){

	bool isBadLivetime=false;
	if(station==2){

		/*
			Anthropogenic movement on Jan 3 2015 between 15:40-15:50 UTC
			Observed in run 4775
		*/
		if( 
			(unixTime>=1389381600 && unixTime<=1389384000) || // from run 2868
			(unixTime>=1420317600 && unixTime<=1420318200) || // from run 4775
			(unixTime>=1449189600 && unixTime<=1449190200) // from run 6507
		)
			{
				isBadLivetime=true;

		}
	}
	else if(station==3){
		// will do something for A3
	}
	return isBadLivetime;

}

/*
	input: station, run number
	output: 0 (is good run), 1 (is bad run)

	function: looks through list of known "bad runs"
			reports if the run you gave it is bad
*/
int isBadRun(int station, int run_num){

	int found;

	/*
	station 2 exclusion
	*/

	vector <double> station2_exclude;

		/*2014*/

			/*
			2014 rooftop pulsing
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			station2_exclude.push_back(3120);
			station2_exclude.push_back(3242);


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

			station2_exclude.push_back(2938); //surface pulsing from m richman (identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14)
			station2_exclude.push_back(2939); //surface pulsing from m richman (identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14)

			/*
			2014 Cal pulser sweep
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			for(int i=3139; i<=3162; i++){ station2_exclude.push_back(i); }
			for(int i=3164; i<=3187; i++){ station2_exclude.push_back(i); }
			for(int i=3289; i<=3312; i++){ station2_exclude.push_back(i); }

			/*
			2014 L2 Scaler Masking Issue
				Cal pulsers sysemtatically do not reconstruct correctly, rate is only 1 Hz
				Excluded because configuration was not "science good"
			*/
			for(int i=3464; i<=3504; i++){ station2_exclude.push_back(i); }

			/*
			2014 Trigger Length Window Sweep
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			for(int i=3578; i<=3598; i++){ station2_exclude.push_back(i); }




		/*2015*/

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
			for(int i=4879; i<=4936; i++){ station2_exclude.push_back(i); }
			for(int i=5210; i<=5277; i++){ station2_exclude.push_back(i); }

			/*
			2015 surface pulsing
				January 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1339 (slide 5)
			*/
			station2_exclude.push_back(4872);
			station2_exclude.push_back(4873);
			station2_exclude.push_back(4876); // Identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14			

			/*
			2015 Pulser Lift
				December 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 2)
				Run number from private communication with John Kelley
			*/
			station2_exclude.push_back(6513);

			/*
			2015 ICL pulsing
				December 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 7)
			*/
			station2_exclude.push_back(6527);


		/*2016*/

			/*
			2016 cal pulser sweep
				January 2015
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2016
			*/
			for(int i=7625; i<=7686; i++){ station2_exclude.push_back(i); }


		/*Other*/

			/*
			D1 Glitches
				Identified by MYL as having glitches after long periods of downtime
			*/
			station2_exclude.push_back(3);
			station2_exclude.push_back(11);
			station2_exclude.push_back(59);
			station2_exclude.push_back(60);
			station2_exclude.push_back(71);			


	/*
	station 3 exclusion
	*/
	
	vector <double> station3_exclude;

		/*2014*/

			/*
			2014 Rooftop Pulser
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			station3_exclude.push_back(2235);
			station3_exclude.push_back(2328);

			/*
			2014 Cal Pulser Sweep
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			for(int i=2251; i<=2274; i++){ station3_exclude.push_back(i); }
			for(int i=2376; i<=2399; i++){ station3_exclude.push_back(i); }

		/*2015*/

			/*
			2015 surface or deep pulsing
				got through cuts
				happened jan 5-6
				waveforms clearly show double pulses or things consistent with surface pulsing
			*/
			station3_exclude.push_back(3811); //deep pulser run
				station3_exclude.push_back(3810); //elminated by proximity to deep pulser run
				station3_exclude.push_back(3820); //elminated by proximity to deep pulser run
				station3_exclude.push_back(3821); //elminated by proximity to deep pulser run
				station3_exclude.push_back(3822); //elminated by proximity to deep pulser run

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
			2015 surface pulsing
				January 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1339 (slide 5)
			*/
			station3_exclude.push_back(3977);
			station3_exclude.push_back(3978);

			/*
			2015 ICL pulsing
				December 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 7)
			*/
			station3_exclude.push_back(6041);

			/*
			Other random
			*/
			station3_exclude.push_back(3977); // looks like ICL events


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
