
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
				if(config==2 || config==4 || config==5){
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
				if(config==2 || config==4){
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
		*/
		if( 
			// Livetime flagged as bad by my me
			(unixTime>=1389381600 && unixTime<=1389384000) || // from run 2868
			(unixTime>=1420317600 && unixTime<=1420318200) || // from run 4775
			(unixTime>=1449189600 && unixTime<=1449190200) || // from run 6507

			//Livetime flagged as bad by my undergrads
			//config 1
			// (unixTime>=1380234000 && unixTime<=1380236400) || // from run 2428 22 hour balloon launch
			// (unixTime>=1382046000 && unixTime<=1382047500) || // from run 2536 22 hour balloon launch
			(unixTime>=1382712900 && unixTime<=1382713500) || // from run 2575 
			(unixTime>=1382972700 && unixTime<=1382973300) || // from run 2589 
			// (unixTime>=1383689100 && unixTime<=1383690900) || // from run 2631 22 hour balloon launch
			(unixTime>=1383884400 && unixTime<=1383886200) || // from run 2642
			(unixTime>=1384060200 && unixTime<=1384061100) || // from run 2652 
			(unixTime>=1384487400 && unixTime<=1384489800) || // from run 2677 
			(unixTime>=1384489980 && unixTime<=1384491060) || // from run 2678 at start may be glitch or continued from 2677
			(unixTime>=1384856520 && unixTime<=1384856640) || // from run 2698 super zoomed in two minute window
			// (unixTime>=1385674200 && unixTime<=1385675100) || // from run 2744 22 hour balloon launch
			(unixTime>=1389381600 && unixTime<=1389383700) || // from run 2868 first of two from run 2868
			(unixTime>=1389398700 && unixTime<=1389400200) || // from run 2868 second of two from run 2868
			(unixTime>=1389665100 && unixTime<=1389666300) || // from run 2884 
			(unixTime>=1393288800 && unixTime<=1393289400) || // from run 3099 
			// (unixTime>=1397856600 && unixTime<=1397858400) || // from run 3442 22 hour balloon launch

			//config 2
			(unixTime>=1376731800 && unixTime<=1376733000) || // from run 2235 

			//conifg 3
			(unixTime>=1400276700 && unixTime<=1400277300) || // from run 3605 mainly looks like glitch at end

			//config 4
			(unixTime>=1409986500 && unixTime<=1409988000) || // from run 4184 
			// (unixTime>=1412026200 && unixTime<=1412027100) || // from run 4301 22 hr balloon
			// (unixTime>=1412285400 && unixTime<=1412288100) || // from run 4316 weird 22hr balloon
			// (unixTime>=1412544600 && unixTime<=1412545500) || // from run 4331 22hr balloon
			// (unixTime>=1412803800 && unixTime<=1412804700) || // from run 4346 22hr balloon
			(unixTime>=1413898200 && unixTime<=1413899100) || // from run 4408
			(unixTime>=1414083900 && unixTime<=1414086000) || // from run 4418 
			(unixTime>=1414350300 && unixTime<=1414351200) || // from run 4434 pt 1
			// (unixTime>=1414358700 && unixTime<=1414359900) || // from run 4434 pt 2 22hr balloon
			(unixTime>=1414674300 && unixTime<=1414674780) || // from run 4452 
			(unixTime>=1414986600 && unixTime<=1414987200) || // from run 4471
			(unixTime>=1415223000 && unixTime<=1415223900) || // from run 4483 
			(unixTime>=1415380500 && unixTime<=1415381400) || // from run 4493 
			(unixTime>=1415558100 && unixTime<=1415559000) || // from run 4503
			(unixTime>=1415742300 && unixTime<=1415743800) || // from run 4513
			(unixTime>=1416207000 && unixTime<=1416212100) || // from run 4541
			(unixTime>=1420978200 && unixTime<=1420978800) || // from run 4814 
			(unixTime>=1416905100 && unixTime<=1416910500) || // from run 4579 two spikes about an hour apart
			// (unixTime>=1416950700 && unixTime<=1416951600) || // from run 4582 22 hour balloon launch
			(unixTime>=1417677000 && unixTime<=1417678200) || // from run 4621  weird and cool
			(unixTime>=1417836000 && unixTime<=1417837500) || // from run 4631
			(unixTime>=1420097100 && unixTime<=1420098300) || // from run 4763 
			(unixTime>=1420293300 && unixTime<=1420294200) || // from run 4774 
			(unixTime>=1420317600 && unixTime<=1420318200) || // from run 4775 
			(unixTime>=1420978200 && unixTime<=1420978800) || // from run 4814 
			(unixTime>=1421024400 && unixTime<=1421025300) || // from run 4817 
			(unixTime>=1421713200 && unixTime<=1421718600) || // from run 4872 looks full of errors and not spiky but could have a spiky
			(unixTime>=1421718000 && unixTime<=1421725800) || // from run 4873 definitely an error but also has spiky boy, part 1 of 2
			(unixTime>=1421733300 && unixTime<=1421733900) || // from run 4873 spiky boy alone but in a run with errors, part 2 of 2
			(unixTime>=1421783400 && unixTime<=1421794200) || // from run 4876 definitely an error but not a spikey boy
			// (unixTime>=1428529800 && unixTime<=1428530700) || // from run 5389 22 hour balloon launch
			(unixTime>=1435623000 && unixTime<=1435623600) || // from run 5801 
			// (unixTime>=1436394000 && unixTime<=1436395200) || // from run 5845 22 hour balloon launch
			(unixTime>=1437601200 && unixTime<=1437602700) || // from run 5915 looks like error at the start
			// (unixTime>=1439933700 && unixTime<=1439934960) || // from run 6048 22 hour balloon launch
			(unixTime>=1440581700 && unixTime<=1440582480) || // from run 6086 
			// (unixTime>=1441489200 && unixTime<=1441490280) || // from run 6137 22 hour balloon launch
			// (unixTime>=1444685400 && unixTime<=1444687080) || // from run 6322 22 hour balloon launch
			// (unixTime>=1445722020 && unixTime<=1445723220) || // from run 6383 22 hour balloon launch
			(unixTime>=1445934900 && unixTime<=1445935500) || // from run 6396 
			(unixTime>=1445960400 && unixTime<=1445961000) || // from run 6397 
			// (unixTime>=1445982120 && unixTime<=1445982900) || // from run 6398 22 hour balloon launch
			(unixTime>=1446165600 && unixTime<=1446166200) || // from run 6408 
			// (unixTime>=1446327300 && unixTime<=1446328200) || // from run 6418 22 hour balloon launch
			(unixTime>=1446607800 && unixTime<=1446608640) || // from run 6433 looks like an error at end
			(unixTime>=1446784200 && unixTime<=1446784800) || // from run 6445 
			// (unixTime>=1476739800 && unixTime<=1476741000) || // from run 8100 22 hour balloon launch
			// (unixTime>=1476999000 && unixTime<=1476999900) || // from run 8114 22 hour balloon launch but barely noticeable
			// (unixTime>=1477258200 && unixTime<=1477259100) || // from run 8129 22 hour balloon launch
			(unixTime>=1477511700 && unixTime<=1477512600) || // from run 8143 weird possible balloon launch
			// (unixTime>=1477950300 && unixTime<=1477951500) || // from run 8168 22 hour balloon launch
			// (unixTime>=1478033400 && unixTime<=1478034000) || // from run 8173 22 hour balloon launch
			// (unixTime>=1478295300 && unixTime<=1478296200) || // from run 8188 22 hour balloon launch
			// (unixTime>=1478728500 && unixTime<=1478729400) || // from run 8213 22 hour balloon launch
			(unixTime>=1479231900 && unixTime<=1479232500) || // from run 8241 

			// config 5
			(unixTime>=1449280500 && unixTime<=1449281100) || // from run 6513
			(unixTime>=1449610200 && unixTime<=1449612000) || // from run 6531
			(unixTime>=1450536000 && unixTime<=1450537200) || // from run 6584
			// (unixTime>=1450906200 && unixTime<=1450907100) || // from run 6606    22hr
			// (unixTime>=1451423700 && unixTime<=1451424600) || // from run 6635   22hr
			(unixTime>=1452008100 && unixTime<=1452009000) || // from run 6669
			// (unixTime>=1452115800 && unixTime<=1452116700) || // from run 6675    22hr
			(unixTime>=1452197700 && unixTime<=1452198600) || // from run 6679
			(unixTime>=1452213600 && unixTime<=1452214200) || // from run 6680
			(unixTime>=1452282000 && unixTime<=1452282600) || // from run 6684 
			(unixTime>=1452298200 && unixTime<=1452298800) || // from run 6685    possible error
			(unixTime>=1452385500 && unixTime<=1452386400) || // from run 6690
			// (unixTime>=1452457800 && unixTime<=1452458700) || // from run 6694   22 hr
			(unixTime>=1452494100 && unixTime<=1452495000) || // from run 6696   possible error
			// (unixTime>=1452545100 && unixTime<=1452545880) || // from run 6700    could be error or 22hr
			// (unixTime>=1452636900 && unixTime<=1452637500) || // from run 6705   could be error or 22hr
			(unixTime>=1452715200 && unixTime<=1452716100) || // from run 6709   possible error
			(unixTime>=1452972300 && unixTime<=1452973440) || // from run 6724   possible error
			// (unixTime>=1453325400 && unixTime<=1453326600) || // from run 6743   22 hr
			(unixTime>=1453408500 && unixTime<=1453409400) || // from run 6747 
			(unixTime>=1453930200 && unixTime<=1453931400) || // from run 6776 
			// (unixTime>=1454535000 && unixTime<=1454536500) || // from run 6818   22 hr
			// (unixTime>=1455746400 && unixTime<=1455747900) || // from run 6889   22 hr
			(unixTime>=1456200900 && unixTime<=1456201800) || // from run 6916
			(unixTime>=1456392600 && unixTime<=1456393800) || // from run 6927
			(unixTime>=1456997400 && unixTime<=1456999200) || // from run 6962
			// (unixTime>=1457559000 && unixTime<=1457560800) || // from run 6994   22 hr
			// (unixTime>=1460842800 && unixTime<=1460844600) || // from run 7119   22 hr
			// (unixTime>=1461620100 && unixTime<=1461621900) || // from run 7161   22 hr
			(unixTime>=1466501400 && unixTime<=1466503200) // from run 7474
			// (unixTime>=1466890200 && unixTime<=1466892000) || // from run 7494   22 hr
			// (unixTime>=1467927600 && unixTime<=1467929700) || // from run 7552   22 hr
			// (unixTime>=1472333400 && unixTime<=1472335200) || // from run 7831   22 hr
			// (unixTime>=1473111300 && unixTime<=1473112800) || // from run 7879    22 hr
			// (unixTime>=1473370500 && unixTime<=1473372900) || // from run 7899   22 hr
			// (unixTime>=1475011500 && unixTime<=1475013600) || // from run 7993   22 hr
			// (unixTime>=1475185200 && unixTime<=1475187900) || // from run 8003 balloon 22hr
			// (unixTime>=1475358000 && unixTime<=1475359800) || // from run 8013 balloon 22hr
			// (unixTime>=1475529900 && unixTime<=1475531400) || // from run 8023 balloon 22hr
			// (unixTime>=1475702700 && unixTime<=1475704200) || // from run 8033 balloon 22hr
			// (unixTime>=1476221400 && unixTime<=1476222300) || // from run 8069 balloon 22hr
			// (unixTime>=1476479700 && unixTime<=1476481800) || // from run 8084 balloon 22hr

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
