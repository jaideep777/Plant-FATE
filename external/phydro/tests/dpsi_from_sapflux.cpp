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
	double psi_s = -1;
	double g;

	phydro::ParPlant P(3e-17, -2, 2);

	phydro::ParEnv E(tc, p, vpd, ppfd/2);
	E.gs_method = phydro::GS_IGF;

	phydro::ParPhotosynth ph(tc, p, kphio, co2, ppfd, fapar, rdark, tc, 25.0);

	auto t1 = std::chrono::high_resolution_clock::now();
	int N = 10;
	for (int i=0; i<N; ++i){
		double dpsi = 0 + i*(10.0)/(N-1);
		double trans = phydro::calc_sapflux(dpsi, psi_s, P, E);
		double dpsi_back = phydro::calc_dpsi_from_sapflux(trans, psi_s, P, E);
		cout << fixed << setprecision(10) << dpsi << "\t" << trans << "\t" << dpsi_back << "\n";

		if (fabs(dpsi - dpsi_back) > 2e-6) return 1;
	}
	auto t2 = std::chrono::high_resolution_clock::now();

	t1 = std::chrono::high_resolution_clock::now();
	N = 50;
	double trans_max = phydro::calc_max_sapflux(psi_s, P, E);
	cout << "Max trans = " << trans_max << '\n';
	for (int i=0; i<N; ++i){
		double trans = trans_max/2 + i*(trans_max)/(N-1);
		double dpsi_back = phydro::calc_dpsi_from_sapflux(trans, psi_s, P, E);
		cout << fixed << setprecision(10) << trans << "\t" << dpsi_back << "\n";

	}
	t2 = std::chrono::high_resolution_clock::now();


	cout << "Time required: " << (std::chrono::duration<double, std::milli> (t2 - t1)).count() << " ms\n";


}


