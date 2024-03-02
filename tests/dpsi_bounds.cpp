#include <iostream>
#include <iomanip>
#include "hyd_transpiration.h"
#include "hyd_analytical_solver.h"
#include <chrono>
using namespace std;

int main(){

	double tc=25;
	double p = phydro::calc_patm(0);
	double vpd = 1000;
	double kphio = 0.087, co2 = 400, ppfd = 1000, fapar = 1, rdark = 0.02;

	double g;

	phydro::ParPlant P(3e-17, -2, 2);

	phydro::ParEnv E(tc, p, vpd, 1000/2);
	E.gs_method = phydro::GS_IGF;

	phydro::ParPhotosynth ph(tc, p, kphio, co2, ppfd, fapar, rdark, tc, 25.0);

	auto t1 = std::chrono::high_resolution_clock::now();
	int N = 10;
	for (int i=0; i<N; ++i){
		double psi_s = -6.0 + i*(6.0)/(N-1);
		phydro::DPsiBounds b = phydro::calc_dpsi_bound(psi_s, P, E, ph, phydro::ParCost(0.1, 1));
		cout << fixed << setprecision(6) << psi_s << "\t" << b.exact << "\t" << b.approx_O2 << "\t" << b.Iabs_bound << "\n";
	}
	auto t2 = std::chrono::high_resolution_clock::now();


	cout << "Time required: " << (std::chrono::duration<double, std::milli> (t2 - t1)).count() << " ms\n";


}


