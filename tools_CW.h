
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

	TGraph *grOut = new TGraph();

	for (int i = 0; i < length1; i++){
		gr1->GetPoint(i, x1, y1);
		gr2->GetPoint(i, x2, y2);
		double phase_diff = y2-y1;

		grOut->SetPoint(i, x1, phase_diff);
	}
	return grOut;

}

double getMedian(TGraph* gr, double lowFreqLimit, double highFreqLimit, double &upper95, double &lower95, double &sigma){
	int npoints = gr->GetN();
	// cout<<"Npoints in get median is "<<npoints<<endl;
	vector < double > vec;
	double x, y;
	for (int i = 0; i < npoints; i++){
		gr->GetPoint(i, x, y);
		// cout << x << " : " << y << endl;
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
		// cout << middle << endl;
		median = (vec[middle]+vec[middle-1])/2.;
	} else {
		middle = int((vec_size-1)/2);
		median = vec[middle];
	}
	int upper95index = middle + int(double(vec_size)/2.*0.95);

	upper95 = vec[upper95index];

	sigma = (upper95 - median)/1.64;
	// cout << "vec_size:upper95index:sigma :::: " << vec_size << " : " << upper95index << " : " << sigma << endl;;

	return median;
}

TGraph *getPhaseVariance( vector<deque<TGraph*> > vdGrPhaseDiff, int runNum, int eventNum, int pol, bool print=false ){
	int numEvents = vdGrPhaseDiff[0].size();
	int numPairs = vdGrPhaseDiff.size();
	// printf("Num pairs is %d \n",numPairs);

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
				// cout << y[i] << " :::: " << real << " : " << im << endl;
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

	if(print){
		char *plotPath(getenv("PLOT_PATH"));
		if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
		TCanvas *c = new TCanvas("","",1100,850);
		gSigmaVarianceAvg->Draw("ALP");
		char save_temp_title[300];
		gSigmaVarianceAvg->SetTitle("Phase Variance Plot");
		gSigmaVarianceAvg->GetYaxis()->SetTitle("Phase variance factor");
		gSigmaVarianceAvg->GetXaxis()->SetTitle("Frequency [MHz]");
		gSigmaVarianceAvg->GetYaxis()->SetRangeUser(-0.5, 3.);
		sprintf(save_temp_title,"%s/phases/Run%d_Ev%d_Pol%d_gSigmaVarianceAvg.png",plotPath,runNum,eventNum,pol);
		c->SaveAs(save_temp_title);
		delete c;
	}

	for ( int i = 0; i < numPairs; i++){
		delete vgPhaseVariance[i];
		delete vgSigmaVariance[i];
	}

	return gSigmaVarianceAvg;

}

vector<double> CWCut_TB(vector <TGraph*> waveforms, vector <TGraph*> baselines, int pol, double dBCut, double dBCutBroad, int station, int num_coinc, vector<int> chan_exclusion_list, int runNum, int eventNum, bool print=false){
	double lowFreqLimit=120.;
	double highFreqLimit=850.;
	double halfrange = (highFreqLimit - lowFreqLimit)/2.;
	double halfway = double(lowFreqLimit + halfrange);
	const int numAnts = 16;
	double deltaTInt = 0.5;

	vector < vector < double> > badFreqs;
	badFreqs.resize(numAnts);
	vector < vector < double> > badFreqsBroad;
	badFreqsBroad.resize(numAnts);

	TGraph *baseline_clone[numAnts];
	vector <TGraph*> newFFTs;
	vector <TGraph*> newBaselines;

	double deltaF_save;

	AraGeomTool *geomTool = AraGeomTool::Instance();

	for(int ant=0; ant<numAnts; ant++){

		double magFFT[2000];
		double magFFTBegin[2000];
		double frequencyArray[2000];
		for(int i=0; i<2000; i++){
			magFFT[i]=0;
			frequencyArray[i]=-1;
		}

		int WaveformLength = 2048; // big, unfortunately...
		// what comes next is a not-so-obvious (imo) way of padding the waveform
		TGraph *chan1Int = FFTtools::getInterpolatedGraph(waveforms[ant],deltaTInt);
		double *getX = chan1Int->GetX();
		double deltaT = getX[1]-getX[0];
		while(chan1Int->GetN() < WaveformLength){
			double lastX = 0.;
			double lastY = 0.;
			chan1Int->GetPoint(chan1Int->GetN()-1,lastX,lastY);
			chan1Int->SetPoint(chan1Int->GetN(), lastX + deltaT, 0 );
		}
		double *getY = chan1Int->GetY();
		int length=chan1Int->GetN();
		FFTWComplex *theFFT=FFTtools::doFFT(length,getY);

		int newLength=(length/2)+1;
		double deltaF=1/(deltaT*length); //Hz
		deltaF*=1e3; //MHz
		deltaF_save=deltaF;

		for(int i=1;i<newLength;i++) {
			if (i==0) frequencyArray[i]=0.;
			if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<=highFreqLimit){
				magFFT[i]+=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;
			}
			else{
				magFFT[i]=-1000;
			}
		}
		delete chan1Int;
		delete [] theFFT;

		for(int i=0;i<newLength;i++) {
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<=highFreqLimit){
				magFFT[i]=10*log10(sqrt(magFFT[i]));
			}
			else magFFT[i]=-1000;
		}

		//need to copy the baselines into new graphs
		double xx, yy;
		double bXreal[newLength];
		double bYreal[newLength];

		baseline_clone[ant] = new TGraph();

		for (int i3 = 0; i3 < newLength; i3++){
			baselines[ant]->GetPoint(i3, xx, yy);
			baseline_clone[ant]->SetPoint(i3, xx, yy);
			bXreal[i3] = xx;
			bYreal[i3] = yy;
		}
		double *bY = baselines[ant]->GetY();
		double *bX = baselines[ant]->GetX();
		int n = baselines[ant]->GetN();

		double *bY_unmodified=baseline_clone[ant]->GetY();
		double *bX_unmodified=baseline_clone[ant]->GetX();
		int n_unmodified=baseline_clone[ant]->GetN();


		// Calculate mean baseline
		//get baseline average so we can bump FFT around
		double mean=0;
		double meanBaseline=0;
		int navg=0;
		int nfirsthalfavg=0;
		int nsecondhalfavg=0;
		double firstHalfMeanBaseline=0;
		double secondHalfMeanBaseline=0;
		double firstHalfMean=0;
		double secondHalfMean=0;

		for (int i=0;i<n;i++){
			if (bX[i]>=lowFreqLimit && bX[i]<highFreqLimit){
				meanBaseline+=bY[i];
				navg++;
			}
			if (bX[i]>=lowFreqLimit && bX[i]<halfway){
				firstHalfMeanBaseline+=bY[i];
				nfirsthalfavg++;
			}
			if (bX[i]>=halfway && bX[i]<highFreqLimit){
				secondHalfMeanBaseline+=bY[i];
				nsecondhalfavg++;
			}
		}
		meanBaseline=meanBaseline/double(navg);
		firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
		secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);

		navg=0;
		nfirsthalfavg=0;
		nsecondhalfavg=0;

		//get average of graph in question
		for (int i=0;i<newLength;i++){
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<highFreqLimit){
				mean+=magFFT[i];
				navg++;
			}
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<halfway){
				firstHalfMean+=magFFT[i];
				nfirsthalfavg++;
			}
			if (frequencyArray[i]>=halfway && frequencyArray[i]<highFreqLimit){
				secondHalfMean+=magFFT[i];
				nsecondhalfavg++;
			}
		}
		mean=mean/double(navg);
		firstHalfMean=firstHalfMean/double(nfirsthalfavg);
		secondHalfMean=secondHalfMean/double(nsecondhalfavg);

		//now bump the average to the baseline average and apply a tilt correction to baseline
		double deltaMean=mean-meanBaseline;
		double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
		double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
		double slope=(deltaMeanFirst-deltaMeanSecond)/halfrange;

		for (int i=0;i<newLength;i++){
			magFFTBegin[i]=magFFT[i];
			magFFT[i]=magFFT[i]-deltaMean;
		}
		for (int ctr=0;ctr<n;ctr++){
			if (bX[ctr]>=lowFreqLimit && bX[ctr]<highFreqLimit){
				bYreal[ctr]=bYreal[ctr]-slope*(bXreal[ctr]-halfway);
			}
		}

		//now see if any peaks are ndB above the baseline.
		double deltaMag[newLength];

		int j;
		for (int i=0;i<newLength;i++){
			if (frequencyArray[i]>lowFreqLimit+2 && frequencyArray[i]<highFreqLimit-2){
				for (j=0;j<n;j++){
					if (bX[j]>frequencyArray[i]) break; // finds bin where frequencies match
				}
				deltaMag[i]=magFFT[i]-bYreal[j]; // changed
			}
			else deltaMag[i]=-1000;
		}

		TGraph *magFFTgraph = new TGraph(2000,frequencyArray,magFFT);
		TGraph *newBaseline = new TGraph(2000,frequencyArray,bYreal);
		newFFTs.push_back(magFFTgraph);
		newBaselines.push_back(newBaseline);

		int index;
		for (int bin = 0; bin < newLength; bin++){
			if (deltaMag[bin] > dBCut){
				badFreqs[ant].push_back(frequencyArray[bin]);
			}
			if (deltaMag[bin] > dBCutBroad){
				badFreqsBroad[ant].push_back(frequencyArray[bin]);
			}
		}
	} //loop over antennas

	vector< vector <double> > freqMatching;
	double freqRangeBroad=40.;
	double freqRangeBroad_bins = freqRangeBroad/deltaF_save; //number of bins spanned by 40 MHz
	vector<double> FreqToNotch;

	for(int i=0; i<numAnts-1; i++){
		int i_pol = geomTool->getStationInfo(station)->getAntennaInfo(i)->polType;
		if(
			i_pol!=pol
			||
			(std::find(chan_exclusion_list.begin(), chan_exclusion_list.end(), i) != chan_exclusion_list.end())
		)
		{
			// check the polarization agreement
			// and to make sure this isn't in the excluded channel list
			// printf("Pol %d: Skipping Antenna %d \n", pol, i);
			continue;
		}
		for(int ii=0; ii<int(badFreqs[i].size()); ii++){ //loop over the bad frequencies for this antenna
			int matchedFreqs=1;
			int matchedFreqsFull=0;
			bool broad_band1=false;
			bool broad_band2=false;
			int broad_freqs1=0;
			int broad_freqs2=0;

			//loop over the other frequencies which were tagged as being above the "broad" threshold in this channel
			for(int ii2=0; ii2<int(badFreqsBroad[i].size()); ii2++){
				//if that frequency is w/in freqRangeBroad window
				if((abs(badFreqs[i][ii] - badFreqsBroad[i][ii2]) < freqRangeBroad)){
					broad_freqs1++;
				}
			}

			//if there are frequencies over the broad threshold w/in 40 MHz
			//then we want to know if the number of contaminated bins/number of bins in 40 MHz < 50%
			//if it's larger, then it's broadband, and we shouldn't touch it!
			if( (freqRangeBroad_bins) != 0){
				if(double(broad_freqs1-1)/double(int(freqRangeBroad_bins)-1) > 0.5){
					broad_band1=true;
				}
				else{
					matchedFreqsFull++;
				}
			}
			else{
				matchedFreqsFull++;
			}

			//now loop over all the other antennas in the array
			for(int j=i+1; j<numAnts; j++){
				int j_pol = geomTool->getStationInfo(station)->getAntennaInfo(j)->polType;
				// cout<<"	Now on antenna j "<<j<<" with pol "<<j_pol<<endl;
				if(
					j_pol!=i_pol
					||
					(std::find(chan_exclusion_list.begin(), chan_exclusion_list.end(), j) != chan_exclusion_list.end())
				)
				{
					// check the polarization agreement
					// and to make sure this isn't in the excluded channel list
					// printf("		One layer deeper. Still on pol %d, now skipping channel %d \n",pol,j);
					continue;
				}
				bool matched_ant = false;
				//check all of their bad frequencies
				for(int jj=0; (jj<int(badFreqs[j].size()) && matched_ant==false); jj++){
					//if their bad frequencies happens within 5 MHz of the first bad frequency
					//we know something is up
					if((abs(badFreqs[i][ii] - badFreqs[j][jj]) < 5.)){

						broad_freqs2=0;
						matchedFreqs++;
						matched_ant=true;
						//now check all of the bad frequencies in this secondary antenna
						for (int jj2 = 0; jj2 < int(badFreqsBroad[j].size()); jj2++){
							//if it's trouble make frequencies are w/in 40 MHz, then we want to record that
							if ((abs(badFreqs[j][jj] - badFreqsBroad[j][jj2]) < freqRangeBroad)){
								broad_freqs2++;
							}
						}
						//if there are frequencies over the broad threshold w/in 40 MHz
						//then we want to know if the number of contaminated bins/number of bins in 40 MHz < 50%
						//if it's larger, then it's broadband, and we shouldn't touch it!
						if(freqRangeBroad_bins!=0){
							if(double(broad_freqs2-1)/double(int(freqRangeBroad_bins)-1) > 0.5){
								broad_band2=true;
							}
						}
						//if the first ant was *not* broadband, and *neither* was this one
						//then we should say "yeah, we found something narrow, please notch me"
						if((broad_band1==false) && (broad_band2==false)){
							matchedFreqsFull++;
						} //were the trouble frequencies for both ants independently not broadband
					} //was the second ant's bad frequency within 5 MHz of the first antennas bad frequency?
				} //loop over second antennas bad freqs
			}//loop over second antenna
			//cout<<"matchedFreqsFull is "<<matchedFreqsFull<<endl;
			//cout<<"Freq "<<badFreqs[i][ii]<<" has matchedFreqsFull of "<<matchedFreqsFull<<endl;
			if(matchedFreqsFull>=num_coinc){
				double new_freq=badFreqs[i][ii];
				for(int k=0; k<FreqToNotch.size(); k++){
					if(abs(new_freq)-FreqToNotch[k] < 0.01)
						new_freq=false;
				}
				if(new_freq)
					FreqToNotch.push_back(badFreqs[i][ii]);
			}
		}//loop over trouble frequencies for antenna 1
	} //loop over antenna 1

	if(print){

		std::vector<TGraph*> newBaselines_6dB;
		for(int i=0; i<16; i++){
			newBaselines_6dB.push_back((TGraph*)newBaselines[i]->Clone());
			for(int samp=0; samp<newBaselines_6dB[i]->GetN(); samp++){
				newBaselines_6dB[i]->GetY()[samp]+=6.;
			}
		}

		char *plotPath(getenv("PLOT_PATH"));
		if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
		TCanvas *c = new TCanvas("","",2*1100,2*850);
		c->Divide(4,4);
		for(int i=0; i<16; i++){
			c->cd(i+1);
			newFFTs[i]->Draw("ALP");
			newBaselines[i]->Draw("same");
			newBaselines[i]->SetLineColor(kRed);
			newBaselines[i]->SetLineWidth(2);
			newFFTs[i]->GetYaxis()->SetRangeUser(15,50);
			newFFTs[i]->GetYaxis()->SetTitle("Power (dB)");
			newFFTs[i]->GetXaxis()->SetTitle("Frequency (MHz)");
			newFFTs[i]->SetLineWidth(2);

			// newBaselines_6dB[i]->Draw("same");
			// newBaselines_6dB[i]->SetLineColor(kBlue);
			// newBaselines_6dB[i]->SetLineWidth(2);

		}
		char save_temp_title[300];
		sprintf(save_temp_title,"%s/trouble_events/Run%d_Ev%d_Pol%d_CWBaseline.png",plotPath,runNum,eventNum,pol);
		c->SaveAs(save_temp_title);
		delete c;
		for(int i=0; i<16; i++) delete newBaselines_6dB[i];
	}

	for(int i=0; i<16; i++){
		delete baseline_clone[i];
		delete newFFTs[i];
		delete newBaselines[i];
	}
	return FreqToNotch;
}

// determine if the baselines themselves have a peak

bool areBaselinesGood(vector<TGraph*> baselines){

	bool areTheseBaselinesGood=true;

	char equation[150];
	sprintf(equation,"([3]*x*x*x + [2]*x*x + [0]*x + [1])");
	char equation_name[16][150];
	TF1 *fit[16];
	double start_of_fit[16];
	double end_of_fit[16];
	for(int chan=0; chan<16; chan++){
		start_of_fit[chan]=150.;
		end_of_fit[chan]=850.;
		sprintf(equation_name[chan],"Fit%d",chan);
		fit[chan] = new TF1(equation_name[chan],equation,start_of_fit[chan],end_of_fit[chan]);
		baselines[chan]->Fit(equation_name[chan],"Q,R");
	}

	int numViolations=0;

	for(int chan=0; chan<16; chan++){
		int numSamps = baselines[chan]->GetN();
		Double_t *theseY = baselines[chan]->GetY();
		Double_t *theseX = baselines[chan]->GetX();
		for(int samp=0; samp<numSamps; samp++){
			double thisX = theseX[samp];
			double thisY = theseY[samp];
			// if(thisX>150. && thisX<850. && (thisX-300.>1.) && (thisX-500.>1.)){
			if(abs(thisX-300.)<2.) continue;
			if(abs(thisX-500.)<2.) continue;
			if(thisX>150. && thisX<800.){
				double expectedY = fit[chan]->Eval(thisX);
				double violation = thisY - expectedY;
				if( violation > 2.){
					// printf("Chan %d, Freq %.2f has %.2f violation \n", chan, thisX, violation);
					numViolations++;
				}
			}
		}
	}

	for(int chan=0; chan<16; chan++) delete fit[chan];

	if(numViolations>0)
		areTheseBaselinesGood=false;

	return areTheseBaselinesGood;

}