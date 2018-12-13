
//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <complex>
#include <deque>

// time
#include <ctime>

#include <AraGeomTool.h>

#include <FFTtools.h>

int getChansfromPair(AraGeomTool * geomTool, int stationNum, int polarization, int pair, int &ch1, int &ch2){
  //  int chan1, chan2;
  int pairnum = 0;
  for (int i = 0; i < 15; i++){
    AraAntPol::AraAntPol_t pol1 = geomTool->getPolByRFChan(i, stationNum);
    if (int(pol1) == polarization){
      for (int i2 = i+1; i2 < 16; i2++){
	AraAntPol::AraAntPol_t pol2 = geomTool->getPolByRFChan(i2, stationNum);
	//	cout << pol1 << " : " << pol2 << endl;
	if (int(pol2) == polarization){
	  if (pairnum == pair){
	    ch1 = i;
	    ch2 = i2;
	    return pairnum;
	  }
	  pairnum++;
	}
      }
    }
  }
  cout << "No chans from pair num" << endl;
  ch1 = -1;
  ch2 = -1;
  return -1;
}



/*
int getPairNum(AraGeomTool *geomTool, int stationNum, int polarization, int polnum1, int polnum2){
 
  
  for (int i = 0; i < 16; i++){
    AraAntPol::AraAntPol_t pol = getPolByRFChan(i, stationNum);
    int polnum = getAntNumByRFChan(i, stationNum);
    if (polnum == polnum1

  }

  }
  
  

}
*/



int getRunNumber_event (const char *runfile){

  string file = string( runfile );
  string chRun = "event";
  size_t foundRun=file.find(chRun);
  string strRunNum = file.substr (foundRun + 5, foundRun + 9);
  int runNum = atoi(strRunNum.c_str());

  return runNum;

}

int padGraph (TGraph *gr, int padLimit){

  int nPointsStart = gr->GetN();

  if (padLimit > nPointsStart){

    double x0, y0;
    gr->GetPoint(nPointsStart-2, x0, y0);
    double x1, y1;
    gr->GetPoint(nPointsStart-1, x1, y1);
    double x_diff = x1-x0;
    
    
    int nNewPoints = padLimit - nPointsStart;
    double x, y;
    for (int i = 0; i < nNewPoints; i++){
      int n_points = gr->GetN();
      gr->GetPoint(n_points-1, x, y);
      gr->SetPoint(n_points, x+x_diff, 0.);
    }
  }

  int nPointsFinal = gr->GetN();
  return nPointsFinal; 
  
}

TGraph * getFFTAmplitude(TGraph *gr, double lowFreqLimit = 0., double highFreqLimit = 1000.);

TGraph * getFFTAmplitude(TGraph *gr, double lowFreqLimit, double highFreqLimit){


  double *getX = gr->GetX();
  double *getY = gr->GetY();

  double deltaT = getX[1] - getX[0];
  int length = gr->GetN();

  FFTWComplex *theFFT = FFTtools::doFFT(length, getY);

  int newLength = (length/2)+1;
  double deltaF = 1/(deltaT*length);
  deltaF*=1e3;

  const int points = newLength;

  double magFFT[points];
  double frequencyArray[points];



  for (int i = 0; i < newLength; i++){
    if (i==0) frequencyArray[i]=0;
    if (i > 0) frequencyArray[i]=frequencyArray[i-1]+deltaF;
    double real2 = theFFT[i].re*theFFT[i].re;
    double im2 = theFFT[i].im*theFFT[i].im;

    if (frequencyArray[i] >= lowFreqLimit && frequencyArray[i] <= highFreqLimit) {
      magFFT[i] = real2+im2;
      magFFT[i] = sqrt(magFFT[i]);
    }
    else {
      magFFT[i] = 0;
    }
  }

  TGraph *grOut = new TGraph(points, frequencyArray, magFFT);

  return grOut;

}

TGraph * getFFTPhase(TGraph *gr, double lowFreqLimit = 0., double highFreqLimit = 1000.);

TGraph * getFFTPhase(TGraph *gr, double lowFreqLimit, double highFreqLimit){


  double *getX = gr->GetX();
  double *getY = gr->GetY();

  double deltaT = getX[1] - getX[0];
  int length = gr->GetN();

  FFTWComplex *theFFT = FFTtools::doFFT(length, getY);

  int newLength = (length/2)+1;
  double deltaF = 1/(deltaT*length);
  deltaF*=1e3;

  const int points = newLength;

  double magFFT[points];
  double frequencyArray[points];
  double phaseFFT[points];

  for (int i = 0; i < newLength; i++){
    if (i==0) frequencyArray[i]=0;
    if (i > 0) frequencyArray[i]=frequencyArray[i-1]+deltaF;
    double real = theFFT[i].re;
    double im = theFFT[i].im;
    double real2 = real*real;
    double im2 = im*im;
    complex<double> fft(real, im);

    if (frequencyArray[i] >= lowFreqLimit && frequencyArray[i] <= highFreqLimit) {
      magFFT[i] = real2+im2;
      magFFT[i] = sqrt(magFFT[i]);
      phaseFFT[i] = arg(fft);
      //      phaseFFT[i] = ata;
    }
    else {
      magFFT[i] = 0;
    }
  }

  TGraph *grOut = new TGraph(points, frequencyArray, phaseFFT);

  delete [] theFFT;

  return grOut;

}

TGraph * getPhaseDifference (TGraph* gr1, TGraph * gr2){
  int length1 = gr1->GetN();
  int length2 = gr2->GetN();

  double x1, y1, x2, y2;

  //  double *frequencyArray = gr1->GetX();
  TGraph *grOut = new TGraph();

  for (int i = 0; i < length1; i++){
    gr1->GetPoint(i, x1, y1);
    gr2->GetPoint(i, x2, y2);
    double phase_diff = y2-y1;

    grOut->SetPoint(i, x1, phase_diff);
  }
  return grOut;

}

//double getMedian(TGraph* gr, double lowFreqLimit = 150., double highFreqLimit = 1000.);
double getMedian(TGraph* gr, double lowFreqLimit, double highFreqLimit, double &upper95, double &lower95, double &sigma){
  int npoints = gr->GetN();
  //cout<<"Npoints in get median is "<<npoints<<endl;
  vector < double > vec;
  double x, y;
  for (int i = 0; i < npoints; i++){
    gr->GetPoint(i, x, y);
    //    cout << x << " : " << y << endl;
    if (x > lowFreqLimit && x < highFreqLimit){
      vec.push_back(y);
    }
  }
  std::sort(vec.begin(), vec.end());
  int vec_size = int(vec.size());
  double median;
  int middle; 
  if (vec_size%2 == 0){
    middle = vec_size/2;
    //cout << middle << endl;
    median = (vec[middle]+vec[middle-1])/2.;
  } else {
    middle = int((vec_size-1)/2);
    median = vec[middle];
  }
  int upper95index = middle + int(double(vec_size)/2.*0.95);

  upper95 = vec[upper95index];

  sigma = (upper95 - median)/1.64;
  //  cout << "vec_size:upper95index:sigma :::: " << vec_size << " : " << upper95index << " : " << sigma << endl;;

  return median;
}

TGraph *getPhaseVariance( vector<deque<TGraph*> > vdGrPhaseDiff ){
  int numEvents = vdGrPhaseDiff[0].size();
  int numPairs = vdGrPhaseDiff.size();

  vector<TGraph*> vgPhaseVariance; vgPhaseVariance.resize(numPairs);
  vector<TGraph*> vgSigmaVariance; vgSigmaVariance.resize(numPairs);
  for (int pairIndex = 0; pairIndex < numPairs; pairIndex++){
    vgPhaseVariance[pairIndex] = new TGraph();
    vgSigmaVariance[pairIndex] = new TGraph();
  }

  for (int pairIndex = 0; pairIndex < numPairs; pairIndex++){
    double x[numEvents];
    double y[numEvents];
    int npoints = vdGrPhaseDiff[pairIndex][0]->GetN();
    
    for (int ii = 0; ii < npoints; ii++){
      complex<double> complex_sum (0,0);
      for (int i = 0; i < numEvents; i++){
	vdGrPhaseDiff[pairIndex][i]->GetPoint(ii, x[i], y[i]);
	double real = cos(y[i]);
	double im = sin(y[i]);
	//			    cout << y[i] << " :::: " << real << " : " << im << endl;
	complex<double> y_com (real, im);
	complex_sum = complex_sum + y_com;
      }
      double phase_variance = 1. - abs(complex_sum)/double(numEvents);
      vgPhaseVariance[pairIndex]->SetPoint(ii, x[0], phase_variance);
    }
    
    double upper95, lower95, sigma;
    double median = getMedian(vgPhaseVariance[pairIndex], 120., 1000., upper95, lower95, sigma);
    int npoints_temp = vgPhaseVariance[pairIndex]->GetN();
    double x0,y0;
    
    
    for(int i = 0; i < npoints_temp; i++){
      vgPhaseVariance[pairIndex]->GetPoint(i,x0,y0);
      if (x0 > 120. && x0 < 1000.){
	int npoints_sigma = vgSigmaVariance[pairIndex]->GetN();
	double sigma_i = (median - y0)/sigma;
	vgSigmaVariance[pairIndex]->SetPoint(npoints_sigma, x0, sigma_i);
      }
    }
  }
  
  TGraph* gSigmaVarianceAvg = new TGraph();
  int npoints_temp = vgSigmaVariance[0]->GetN();
  for (int i = 0; i < npoints_temp; i++){
    double average = 0.;
    double x1,y1;
    for (int pairIndex = 0; pairIndex < numPairs; pairIndex++){
      vgSigmaVariance[pairIndex]->GetPoint(i,x1,y1);
      average = average + y1;
    }
    average = average / (double)numPairs;
    gSigmaVarianceAvg->SetPoint(i, x1, average);
  }
  
  for ( int i = 0; i < numPairs; i++){
    delete vgPhaseVariance[i];
    delete vgSigmaVariance[i];
  }
  
  return gSigmaVarianceAvg;

}
