
#ifndef FILTEREVENT_H
#define FILTEREVENT_H


//class TGraph;
//class TH2D;

///#include <vector>;


class filterEvent{
double qualityParameter;
std::vector< std::vector < double > > hits;
double chNoiseEnergy[16];
double chNoiseEnergyRMS[16];
double chSignalEnergy[16];
double chSignalMaxEnergyPosition[16];
double ch_noise_std_dev[16];
double ch_noise_average[16];
std::vector<std::vector< double > > energyEnvelope;
std::vector<std::vector< double > > energyTime;
int totalHitCount;
int binnedHitCount;
double blockCount;
public:
filterEvent();
double getHorizontalHitDistance(int hit1, int hit2);
double getQualityParameter(std::vector<TGraph*> gr_in, std::vector<std::vector<double> > ant_loc, int stationId, double * qualArray);
double gen_trig_histo_single(std::vector<TGraph *> gr, TH2D *trig_pat, double *ch_std_dev, double *ch_average, double th_factor, int sum_time, int stationId);

int gen_noise_std_dev(std::vector<TGraph *> gr, double *ch_noise_std_dev, double *ch_noise_average, int sum_time);

double pattern_check(TH2D *trig_pat, double *pat_check_result, int cut_value, std::vector< std::vector< double > > antloc);

int inter_string_check(TH2D *trig_pat, int *same_speed_t, int *same_speed_d1, int *same_speed_d2, int *same_speed_v, int *same_speed_h, int *c_s, int *c_d, std::vector< std::vector< double > > antloc);

double getMaxEnergy(std::vector<TGraph *> gr);

double getPolarization(int stationId);
TGraph * getEnergyGraph(int stationId, int channel);
double getAverageNoiseEnergy(int stationId, int channel);
double getAverageNoiseEnergyRMS(int stationId, int channel);
double getChannelMaxEnergy(int stationId, int channel);
double getSNR(int stationId);

double getTotalHitCountValue();
double getBinnedHitCountValue();

};


#endif
