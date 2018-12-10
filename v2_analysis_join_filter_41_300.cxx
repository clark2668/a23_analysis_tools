////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_join_filter_41_300.cxx 
////	A23 diffuse, merge filter, 41m, and 300m files to joined files
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////


//Includes
#include <iostream>
#include <string>
#include <sstream>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

#include "inputParameters.h"

using namespace std;

int main(int argc, char **argv)
{

	int recoBinCalpulser = 6;
	int recoBinSelect = 19;

	stringstream ss;

	if(argc<2) {
		std::cout << "Usage\n" << argv[0] << " <base file name>";
		std::cout << "e.g.\n" << argv[0] << " ./outputs http://www.hep.ucl.ac.uk/uhen/ara/monitor/root/run1841/event1841.root default_pedestal.txt\n";
		return 0;
	}

	ss.str("");
	ss << argv[1] << "_joined_bins_" << recoBinCalpulser << "_" << recoBinSelect << ".root";

	TFile *OutputFile = TFile::Open(ss.str().c_str(), "RECREATE");

	TFile *fOpen_filter;
	TTree *inputSettingsTree_filter; 
	TTree *inputTree_filter; 

	ss.str("");
	ss << argv[1] << "_filter.root";
	fOpen_filter = TFile::Open(ss.str().c_str());
	if(!fOpen_filter) {
		std::cerr << "Can't open file : filter" << "\n";
	return -1;
	}

	inputTree_filter = (TTree*) fOpen_filter->Get("OutputTree");
	if(!inputTree_filter) {
		std::cerr << "Can't find OutputTree: filter"  << "\n";
		return -1;
	}

	int nEvents = inputTree_filter->GetEntries();

	OutputFile->cd();  
	ss.str("");
	ss << "OutputTree_filter";
	inputTree_filter->CloneTree()->Write(ss.str().c_str());

	TFile *fOpen_reco[35];
	TTree *inputSettingsTree_reco[35];
	TTree *inputTree_reco[35];
  
	for (int i = 0; i < 35; i++){
		if (i == recoBinSelect || i == recoBinCalpulser){

			ss.str("");
			ss << argv[1] << "_recoRadius_" << radii[i] << ".root";
			fOpen_reco[i] = TFile::Open(ss.str().c_str());
			if(!fOpen_reco[i]) {
				std::cerr << "Can't open file : "<< i << "\n";
				return -1;
			}
      
			inputSettingsTree_reco[i] = (TTree*) fOpen_reco[i]->Get("OutputSettingsTree");
			if(!inputSettingsTree_reco[i]) {
				std::cerr << "Can't find OutputSettingsTree: " << i << "\n";
				return -1;
			}
      
			inputTree_reco[i] = (TTree*) fOpen_reco[i]->Get("OutputTree");
			if(!inputTree_reco[i]) {
				std::cerr << "Can't find OutputTree: " << i << "\n";
				return -1;
			}
      
			int nEvents_reco = inputTree_reco[i]->GetEntries();
			if (nEvents != nEvents_reco){
				cerr << "Event numbers don't match: RecoBin " << i << endl;
				return -1;
			}
      
			OutputFile->cd();
			if (i==recoBinSelect){
				inputSettingsTree_reco[i]->CloneTree()->Write();
			}
      
			ss.str("");
			ss << "OutputTree_recoRadius_" << i;
			inputTree_reco[i]->CloneTree()->Write(ss.str().c_str());

			fOpen_reco[i]->Close();
		}
	}
	OutputFile->Close();
	fOpen_filter->Close();
}
