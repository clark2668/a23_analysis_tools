//////////////////////////////////////////////////////////////////////////////
/////  tools_RecoFns.h                                                   /////
/////                                                                    /////
/////  Description:                                                      /////
/////     Functions extracting info from interferometric maps            /////
/////  Author: CGP                                                       /////
//////////////////////////////////////////////////////////////////////////////

//includes
#include "TH2D.h"

using namespace std;

void getCorrMapPeak_wStats( TH2D *theCorrMap_input, int &peakTheta, int &peakPhi, double &peakCorr, double &minCorr, double &meanCorr, double &rmsCorr, double &peakSigma) {
	
	int peakZ;
	
	theCorrMap_input->GetMaximumBin(peakPhi, peakTheta, peakZ);
	peakCorr = theCorrMap_input->GetMaximum();
	
	peakTheta = peakTheta - 90;
	peakPhi = peakPhi - 180;

	minCorr = theCorrMap_input->GetMinimum();
	meanCorr = theCorrMap_input->GetMean(3);
	rmsCorr = theCorrMap_input->GetRMS(3);

	int nCells = theCorrMap_input->GetSize() - theCorrMap_input->GetNbinsX()*2 - theCorrMap_input->GetNbinsY()*2 - 4;

	double stats[6];
	theCorrMap_input->GetStats(stats);
	meanCorr = stats[0]/(double)nCells;
	
	rmsCorr=0;
	peakSigma=0;

	double sumOfDeviation2 = 0.;
	for (int theta=1; theta<=180; theta++) {
		for (int phi=1; phi<=360; phi++) {
			double corrValue = theCorrMap_input->GetBinContent(phi,theta);
			double deviation = corrValue-meanCorr;
			sumOfDeviation2 = sumOfDeviation2 + pow(deviation,2.);
		}
	}
	rmsCorr = sqrt(sumOfDeviation2/(double)nCells);
	peakSigma = (peakCorr-meanCorr)/rmsCorr;

}


void getCorrMapPeak( TH2D *theCorrMap_input, int &peakTheta, int &peakPhi, double &peakCorr) {
	
	int peakZ;
	
	theCorrMap_input->GetMaximumBin(peakPhi, peakTheta, peakZ);
	peakCorr = theCorrMap_input->GetMaximum();
	
	peakTheta = peakTheta - 90;
	peakPhi = peakPhi - 180;
}

/* deprecated
void getCorrMapPeak( TH2D *theCorrMap_input, int &peakTheta, int &peakPhi, double &peakCorr) {
  double corrPeak = -1.;
	for (int theta=1; theta<=180; theta++) {
		for (int phi=1; phi<=360; phi++) {
			double corrValue = theCorrMap_input->GetBinContent(phi,theta);
			if ( corrValue > corrPeak ) {
				corrPeak = corrValue;
				peakTheta = theta - 90;
				peakPhi = phi - 180;
				peakCorr = corrPeak;
			}
		}
	}
}
*/