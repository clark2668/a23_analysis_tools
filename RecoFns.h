////////////////////////////////////////////////////////////////////////////// 
/////  WaveformFns.h                                                     /////  
/////                                                                    /////  
/////  Description:                                                      ///// 
/////     Functions for making a variety of cuts on Ara Data             ///// 
/////  Author: CGP, adapted from ACG, class structure from RJN           ///// 
////////////////////////////////////////////////////////////////////////////// 
//System includes                                                                                                               
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <sstream>
#include <functional>
#include <numeric>
#include <deque>
using namespace std;

//AraRoot Includes                                                                                                              
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"
#include "AraAntennaInfo.h"
//#include "RayTraceCorrelator_test.h"

#include "FFTtools.h"

//ROOT Includes                                                                                                                 
#include "TROOT.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TButton.h"
#include "TGroupButton.h"
#include <TGClient.h>
#include "TStyle.h"
#include "TPostScript.h"
#include "TTree.h"
#include "math.h"
#include "TText.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "Math/Interpolator.h"
#include "TImage.h"
#include "TMarker.h"
#include "TStyle.h"

#include "Constants.h"

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
