
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

using namespace std;

// compute the log-likelihood 
double ReturnLogL_highN( double N_bin, double U_exp ) {
    double LogL; // -2 ln( L )
	if ( N_bin > 0. ) {
		LogL = -2.* ( N_bin * log( U_exp / N_bin ) - U_exp + N_bin );
	}
	else {
		LogL = -2.* ( -U_exp );
	}
	return LogL;
}
