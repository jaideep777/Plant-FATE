#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "incgamma.h"
#include <gsl/gsl_sf_gamma.h>
// #include "hyd_transpiration.h"

int main(){
	for (double a : vector<double>({0.0001, 0.5, 1, 1.5, 2, 2.5})){
		for (double x=0; x<10; x=x+0.5){
			double res_gsl = gsl_sf_gamma_inc(a, x);
			double res_inc = gammainc(a,x);
			cout << "Gamma(" << a << ", " << x << ") = " << res_gsl << " | " << res_inc << "\n";
			if (fabs(res_gsl-res_inc)>1e-6) return 1;
		}
	}

	// for (double b : vector<double>({0.5, 1, 1.5, 2, 2.5})){
	// 	for (double dpsi=0.1; dpsi<2; dpsi=dpsi+0.1){
	// 		double res = phydro::integral_P_analytical(dpsi, -1, -1.5, b);
	// 		cout << "T(" << b << ", " << dpsi << ") = " << res << "\n";
	// 	}
	// }
	
	return 0;
}

