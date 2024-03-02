#include <iostream>
#include "hyd_transpiration.h"
#include "temperature_dependencies_physical.h"
#include <chrono>
using namespace std;

int main(){

	double tc=25;
	double p = phydro::calc_patm(0);
	double vpd = 1000;
	
	double g;
	phydro::ParPlant P(3e-17, -2, 2);

	phydro::ParEnv E(tc, p, vpd, 1000/2);
	E.gs_method = phydro::GS_IGF;

	auto t1 = std::chrono::high_resolution_clock::now();
	int N = 10;
	for (int i=0; i<N; ++i){
		double psi_s = -6.0 + i*(6.0)/(N-1);
		double Q = phydro::calc_sapflux(1, psi_s, P, E);
		g = phydro::calc_gs_from_Q(Q, psi_s, P, E);
		cout << psi_s << "\t" << g << "\n";
	}
	auto t2 = std::chrono::high_resolution_clock::now();

	cout << g << endl;
	cout << "Time required: " << (std::chrono::duration<double, std::milli> (t2 - t1)).count() << " ms\n";


	if (abs(g-0.1116382) < 1e-5) return 0;
	else return 1;
}


