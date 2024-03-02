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

	phydro::ParPhotosynth ph_test(28, p, kphio, co2, ppfd, fapar, rdark, 19, 25.0);
	ph_test.print();

	double q1 = phydro::calc_sapflux(1, -1, P, E);
	double g1 = phydro::calc_gs_from_Q(q1, -1, P, E);
    double g1p = phydro::calc_gsprime(1, g1, -1, P, E);
    double x = calc_x_from_dpsi(1, g1p, ph, phydro::ParCost(0.1, 1));
    double J = calc_J(g1, x, ph)-4*ph.phi0*ph.Iabs;
	cout << "dFdx f2: " << q1 << " " << g1 << " " << g1p << " " << x << " " << J << "\n";

	auto t1 = std::chrono::high_resolution_clock::now();
	int N = 10;
	for (int i=0; i<N; ++i){
		double psi_s = -6.0 + i*(6.0)/(N-1);
		phydro::DPsiBounds b1 = phydro::calc_dpsi_bound(psi_s, P, E, ph, phydro::ParCost(0.1, 1));
		phydro::DFDX b = phydro::dFdx(b1.Iabs_bound*0.5, psi_s, P, E, ph, phydro::ParCost(0.1, 1));
		cout << fixed << setprecision(6) << psi_s << "\t" << b1.exact << "\t" << b1.Iabs_bound << "\t" << b.dJ_dchi << "\t" << b.djmax_dJ << "\t" << b.dPdx << "\t" << b.J << "\n";
	}
	auto t2 = std::chrono::high_resolution_clock::now();


	cout << "Time required: " << (std::chrono::duration<double, std::milli> (t2 - t1)).count() << " ms\n";


}


