#include <iostream>
#include "hyd_photosynthesis.h"
using namespace std;


int check(double x, double ref, double err=1e-5){
	if (abs(x-ref) < err) return 0;
	else return 1;
}


int main(){

	double kphio = 0.05;        // quantum yield efficiency
	// c_molmass <- 12.0107 # molar mass, g / mol

	// Define environmental conditions
	double tc = 25;             // temperature, deg C
	double ppfd = 400;          // umol/m2/s
	double vpd  = 1000;         // Pa
	double co2  = 400;          // ppm
	double elv  = 0;            // m.a.s.l.
	double fapar = 0.7;         // fractioni
	double rdark = 0;

	double gs = 0.1;
	
	phydro::ParPhotosynth par_photosynth(tc, phydro::calc_patm(elv), kphio, co2, ppfd, fapar, rdark, tc, 25.0);

	auto ac = phydro::calc_assim_rubisco_limited(gs, 30, par_photosynth);
	auto aj = phydro::calc_assim_light_limited(gs, 100, par_photosynth);

	auto al = phydro::calc_assimilation_limiting(30, 100, gs, par_photosynth);

	cout << ac.a << " " << ac.ci << " " << aj.a << " " << aj.ci << " " << al.a << " " << al.ci << endl;

	int nerr = 0;
	nerr += check(ac.a,  8.1325907834);
	nerr += check(ac.ci, 32.289652);
	nerr += check(aj.a,  6.2731149469);
	nerr += check(aj.ci, 34.173766);
	nerr += check(al.a,  6.2731149469);
	nerr += check(al.ci, 34.173766);

	return nerr;
}


