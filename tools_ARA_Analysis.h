////////////////////////////////////////////////////////////////////////////// 
/////  WaveformFns.h                                                     /////  
/////                                                                    /////  
/////  Description:                                                      ///// 
/////     Functions for making a variety of cuts on Ara Data             ///// 
/////  Author: CGP, adapted from ACG, class structure from RJN           ///// 
////////////////////////////////////////////////////////////////////////////// 


#ifndef ARA_ANALYSIS_H
#define ARA_ANALYSIS_H

class ARA_Analysis{

  vector<vector<vector<vector<int > > > > Pairs;
  vector<vector<vector<vector<int > > > > Faces;
  vector<vector<int > > Pairs_for_corr_V;
  vector<vector<int > > Pairs_for_corr_H;

 public:  
  ARA_Analysis();
  ARA_Analysis(int stationID);

  vector<vector<vector<vector<int> > > > setupPairs();

  double getRMS( TGraph *plot, int numPointsToInclude);
  double getRMS( double *array, int bin);
  void getAbsMaximum(TGraph *plot, double &x_max, double &y_max);
  void getAbsMaximum(vector<TGraph*> graphs, vector<double> &xs, vector<double> &ys );
  void getMaximum(TGraph * plot, double &x_max, double &y_max);

  double integrateBinPower( TGraph *plot, int numBinsToIntegrate, vector<double> &integratedBins);  
  vector<TGraph*> makeIntegratedBinPowerGraphs(vector<TGraph*> graphsInput, int numBinsToIntegrate, string xlabel, string ylabel, vector<string> titles);  


  double getWaveformVariance(double *array, int bin);
  double getWaveformVariance( TGraph *plot, int numPointsToInclude);

  double integrateAllPower(TGraph *plot);

  double getAveragePower(TGraph *plot);

  vector<int> getDividers(TGraph *plot, int divisions);

  double integrate(TGraph *plot, int firstBin, int lastBin);

  double getAverage(TGraph *plot, int firstBin, int lastBin);

  double getDeviation(TGraph *plot, int firstBin, int lastBin, double average);

  void getAverageDeviation_Divided(TGraph *plot, int divisions, double &average_out, double &deviation_out);

  void getVAveragesDeviations_Divided(vector<TGraph*> vPlots, int divisions, vector<double> &vAverages, vector<double> &vDeviations);

  double integratePowerFirstNBins(TGraph *plot, int nBins);

  double getAveragePowerFirstNBins(TGraph *plot, int nBins);

  vector<TGraph*> makeAvgVoltageFromIntPower(vector<TGraph*> graphsInput, int numBinsToIntegrate, string xlabel, string ylabel, vector<string> titles);

  void getRMS(vector<TGraph*> graphs, vector<double> &vRMS, int numPointsToInclude);

  void getWaveformVarianceSamples(vector<TGraph*> graphs, int numPointsToInclude, vector<double> &vVariance, vector<double> &vSamples);

  void addVarianceSamples(vector<double> vWaveformVariance, vector<int> vSamples, vector<double> vWaveformVarianceTotal, vector<long int> vSamplesTotal);

  vector<double> getTotalPowers(vector<TGraph*> graphs);

  vector<double> getAveragePowers(vector<TGraph*> graphs);

  vector<double> getAveragePowersFirstNBins(vector<TGraph*> graphs, int nBins);

  vector<vector<vector<vector<int> > > > setupPairs(int stationId);

  vector<vector<vector<double> > > getHitTimeDifferences(vector<double> hitTimes, vector<double> peakIntPowers, vector<vector<vector<vector<int> > > > Pairs);
  

  vector<vector<vector<double> > > getHitTimeWeight_LesserPower(vector<double> hitTimes, vector<double> peakIntPowers, vector<vector<vector<vector<int> > > > Pairs);
  
  vector<vector<vector<double> > > getHitTimeWeight_LesserFractionalPower(vector<double> hitTimes, vector<double> peakIntPowers, vector<double> totalPowers, vector<vector<vector<vector<int> > > > Pairs);

  vector<vector<vector<double> > > getHitTimeWeight_LesserRelativePower(vector<double> hitTimes, vector<double> peakIntPowers, vector<double> averagePowers, vector<vector<vector<vector<int> > > > Pairs);
  

  vector<vector<vector<double> > > getHitTimeWeight_LesserRelativeRMS(vector<double> hitTimes, vector<double> peakIntRMS, vector<double> waveformRMS, vector<vector<vector<vector<int> > > > Pairs);
  
  vector<TH1D*> makeHitTimeDist(vector<vector<vector<double> > > delays, int polarization);

  double getMean(vector<double> v);
		      
  double getRms(vector<double> v);

  vector<vector<double> > getRmsDelays(vector<vector<vector<double> > > delays);

  double getRmsPol(vector<vector<vector<double> > > delays, int polarization);

  double getWeightedSum(vector<double> v, vector<double> weights);

  double getWeightedMean(vector<double> v, vector<double> weights);

  double getWeightedVariance(vector<double> v, vector<double> weights);

  double getWeightedRms(vector<double> v, vector<double> weights);
  
  double getWeightedRmsPol(vector<vector<vector<double> > > delays, vector<vector<vector<double> > > weights, int polarization);
  
  vector<vector<vector<double> > > getVVVToNPower(vector<vector<vector<double> > > vvvIn, double power);
  
  vector<vector<vector<double> > > getHitTimeMagnitude_Lesser(vector<double> hitTimes, vector<double> magnitude, vector<vector<vector<vector<int> > > > Pairs);
  
  bool isUsablePattern(vector<double> magnitudes, double threshold);
  
  double getWeightedSum_Threshold(vector<double> v, vector<double> weights, vector<double> magnitudes, double threshold);

  double accumulateWeights_Threshold(vector<double> weights, vector<double> magnitudes, double threshold);

  double getWeightedMean_Threshold(vector<double> v, vector<double> weights, vector<double> magnitudes, double threshold);

  void getWeightedVariance_Threshold(vector<double> v, vector<double> weights, vector<double> magnitudes, double threshold, double &sumOfSquares_out, double &sumOfWeights_out);

  double getWeightedRmsPol_Threshold(int polarization, vector<vector<vector<double> > > delays, vector<vector<vector<double> > > weights, vector<vector<vector<double> > > magnitudes,  double threshold );

  double getLeastDifferenceFromTetrad(vector<vector<vector<double> > > delays);

  int getHitsUnderTimeFromTetrad(int polarization, double timeLimit, vector<vector<vector<double> > > delays);

  int getHitsUnderTimeFromTetrad(int polarization, double timeLimit, vector<vector<vector<double> > > delays, vector<vector<vector<double> > > weights);

  vector<vector<vector<double> > > getHitTimeDifferencesFromTetrad(vector<double> hitTimes, vector<double> peakIntPowers, vector<vector<vector<vector<int> > > > Pairs);


  vector<double> getVSigmas(vector<double> peaks, vector<double> averages, vector<double> variances);

  vector<vector<vector<double> > > getVVVSigmas(vector<double> vSigmas, vector<vector<vector<vector<int> > > > Pairs);

  double getLesser(double input1, double input2);

  void getMinimum(int length, double *array, int &index, double &minimum);

  double* getCorrelation_NoNorm(int length, double *oldY1, double *oldY2);

  TGraph *getCorrelationGraph_WFweight(TGraph *gr1, TGraph *gr2);

  vector<TGraph*> getCorrelationGraphs_wBestTimes(vector<TGraph*> grIn, vector<vector< int > > pairs, vector<double> &bestTimes, vector<double> &bestCorrs);

  vector<TGraph*> getCorrelationGraphs(vector<TGraph*> grIn, vector<vector< int > > pairs );

  vector<TGraph*> getHilbertGraphs_wBestTimes(vector<TGraph*> grIn, vector<double> &bestTimes, vector<double> &bestCorrs);

  vector<TGraph*> getHilbertCorrelationGraphs_wBestTimes(vector<TGraph*> grIn, vector<vector< int > > pairs, vector<double> &bestTimes, vector<double> &bestCorrs);

  void setupCorrelationPairs(int StationID, vector<vector<int> > &pairs_V, vector<vector<int> > &pairs_H);


  vector<double> getRms_Corr(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc);

  vector<double> getRms_Corr_Threshold(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc, vector<double> BestCorrs, double CorrThreshold);

  double  getPolarizationRatio(vector<TGraph*> waveforms, vector<int> polarizations);

  void getThirdVPeakOverRMS(vector<double> vPeakOverRms, vector<int> polarizations, vector<double> &ThirdVpeakOverRms);

  double getRms_Remove1(vector<double> v);

  double getRms_Remove1(vector<double> v, int &dropped_pair);

  vector<double> getRms_Corr_Regular_Threshold_Remove1(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc, vector<double> BestCorrs, double CorrThreshold);

  vector<double> getRms_Corr_Regular_Threshold_Remove1(vector<double> BestTimes, vector<vector<int> > pairs_from_correlation, vector<vector<vector<vector<int> > > > pairs_from_faces, vector<vector<double> > ant_loc, vector<double> BestCorrs, double CorrThreshold, vector<int> &dropped_pairs);

  vector<vector<vector<vector<int> > > > setupFaces(int stationId);

vector<double> getRms_Faces(
			    vector<double> BestTimes,
			    int polarization,
			    vector<vector<int> > pairs_from_correlation, 
			    vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, 
			    vector<vector<double> > ant_loc);

vector<double> getRms_Faces(
			    vector<double> BestTimes,
			    // vector<double> BestSignal,
			    int polarization,
			    vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, 
			    vector<vector<double> > ant_loc);

 bool isCutOnFaceError(int nFaces, double* faceCuts, vector<double> faceErrors);

vector<double> getRms_Faces_Thresh(
				   vector<double> BestTimes, // vector of arrival times at channels
				   vector<double> BestSignals, // vector of signal strengths to compare with threshold found at arrival times
				   double threshold, // threshold of signal strengths needed to be considered of sufficient strength to be included in face
				   int polarization, // polarization of antennas to check
				   vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, // vector of antenna pairs for faces
				   vector<vector<double> > ant_loc // vector of antenna locations 
				   );

 void getAbsMaximum_N(TGraph* plot, int nPeaks, double timeApart, vector<double> &xs, vector<double> &ys);

 void getAbsMaximum_N(vector<TGraph*> graphs, int nPeaks, double timeApart, vector<vector<double> > &xs, vector<vector<double> > &ys);

vector<double> getRms_Faces_Thresh_N(
				   vector<vector<double> > BestTimes, // vector of arrival times at channels
				   vector<vector<double> > BestSignals, // vector of signal strengths to compare with threshold found at arrival times
				   double threshold, // threshold of signal strengths needed to be considered of sufficient strength to be included in face
				   int polarization, // polarization of antennas to check
				   vector<vector<vector<vector<int> > > > pairs_pol_faces_pair_ind, // vector of antenna pairs for faces
				   vector<vector<double> > ant_loc // vector of antenna locations 
				     );

 double getPhase(FFTWComplex &theNum);

 TGraph *makeRawPhase(TGraph *grWave);

 TGraph *getGraphDifference(TGraph *gr1, TGraph* gr2);

 TGraph *getGraphSum(TGraph *gr1, TGraph* gr2);

 TGraph *makeSpectrum_mVPerRootHz(TGraph *grWave);
 
 void getWavefrontRMS(vector<TGraph*> graphs_raw, vector<double> waveformRMS, 
		       int stationID, vector< vector<double> > ant_loc,
		       double threshold,
		       double &wavefrontRMS_V, double &wavefrontRMS_H,
		       double interpolationTimeStep,
		       double intergrationTime,
		       int numSearchPeaks,
		       double peakSeparation
		      );

 void getWavefrontRMSVector(vector<TGraph*> graphs_raw, vector<double> waveformRMS, 
		       int stationID, vector< vector<double> > ant_loc,
		       vector<double> threshold,
		       vector<double> &vWavefrontRMS_V, vector<double> &vWavefrontRMS_H,
		       double interpolationTimeStep,
		       double intergrationTime,
		       int numSearchPeaks,
		       double peakSeparation
		      );


 void saveGraph(TGraph* gr1, string filename);
 void save2Graphs(TGraph* gr1, TGraph* gr2, string filename);
 void save16Graphs(vector<TGraph*> grs, string filename);
 void saveNGraphs(vector<TGraph*> grs, string filename, int numGraphs);
 void saveNGraphs_2Tlines(vector<TGraph*> grs, string filename, int numGraphs, vector<double> averages, vector<double> deviations);
 
 void saveHist(TH1* hist, string filename);
 void saveHistWithFit(TH1* hist, TF1* myFit, string filename);
 void saveHistOverlap(vector<TH1D*> hists, string filename);
 void saveHistOverlapPolarized(vector<TH1D*> histsV, vector<TH1D*> histsH, string filename, int logOpt=0);
 void saveHist2D(TH2* hist, string filename);
 void saveHist2D(TH2* hist, string filename, int logOpt);
 void saveNHists2D(vector<TH2D*> hists, string filename_base);

 void setHistLineColors(vector<TH1D*> hists);
 
 void deleteGraphVector(vector<TGraph*> graphs);
 void deleteHistVector(vector<TH1D*> hists);

 vector<TGraph*> makeGraphsFromRF(UsefulAtriStationEvent* realAtriEvPtr, int numGraphs, string xlabel, string ylabel, vector<string> titles);
 vector<TGraph*> makeInterpolatedGraphs(vector<TGraph*> graphsIn, double intTimestep, string xlabel, string ylabel, vector<string> titles);
 vector<TGraph*> makePaddedGraphs(vector<TGraph*> graphsIn, bool isSoftTrigger, string xlabel, string ylabel, vector<string> titles);
 vector<TGraph*> makePowerSpectrumGraphs(vector<TGraph*> graphsIn, string xlabel, string ylabel, vector<string> titles);
 vector<TGraph*> makeVoltageSpectrumGraphs(vector<TGraph*> graphsIn, string xlabel, string ylabel, vector<string> titles);
 vector<TGraph*> makePhaseGraphs(vector<TGraph*> graphsIn, string xlabel, string ylabel, vector<string> titles);

 vector<TGraph*> initializeGraphVector(int size);
 vector<TH1D*> initializeHistVector(int size, string name, int bins, double lowBound, double upperBound);
 vector<TH2D*> initializeHistVector2D(int size, string name, int bins_x, double lowBound_x, double upperBound_x, int bins_y, double lowBound_y, double upperBound_y);

 int getRunNum(char* runfile);
 int getrunNum(char* runfile);

 string getProcessedFilename(int stationID, char* outputdir, char* runfile );
 string getProcessedFilename_filter(int stationID, char* outputdir, char* runfile );
 string getProcessedFilename_recoRadius(int stationID, char* outputdir, char* runfile, int radius);
 string getNoiseFilename(int stationID, char* outputdir, char* runfile );
 string getRunSummaryFilename(int stationID, char* outputdir, char* runfile );

 vector<TGraph*> makeDifferenceGraphs(vector<TGraph*> graphsIn, int graphBegin, int graphEnd, string xlabel, string ylabel, vector<string> titles);
 vector<TGraph*> makeSumGraphs(vector<TGraph*> graphsIn, int graphBegin, int graphEnd, string xlabel, string ylabel, vector<string> titles);


};

#endif
