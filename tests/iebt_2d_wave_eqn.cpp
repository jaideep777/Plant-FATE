#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
#include "individual_base.h"
using namespace std;

#include "wave_model_2d.h"

std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

int main(){

	Species<Wave> spp;
	WaveEnv E;

	Solver S(SOLVER_IEBT);
	S.control.ode_ifmu_stepsize = 0.02;
	S.setEnvironment(&E);
	S.addSpecies({100, 100}, {0, 0}, {10,10}, {false, false}, &spp, 0, -1);
	S.print();

	ofstream fout1;
	fout1.open("iebt2d_u.txt");
	S.species_vec[0]->printCohortVector(fout1);
	fout1.close();

	// for (int i=0; i<50; ++i){
	// 	S.copyCohortsToState();
	// 	S.stepU_iFMU(0, S.state, S.rates, 0.02);
	// 	S.copyStateToCohorts(S.state.begin());
	// 	// S.print();
	// }

	ofstream foutx("iebt_wave_x.txt");
	ofstream fouty("iebt_wave_y.txt");

	ofstream fout("iebt_wavemodel_equil.txt");

	for (double t=0; t <= 1; t=t+0.1) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << " " << S.u0_out(t)[0] << "\n";
		
		vector<double> breaks = myseq(0,10,201);
		vector<double> vx = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		for (auto x : vx) foutx << x << "\t";
		foutx << endl;

		vector<double> vy = S.getDensitySpecies1D(0, 1, breaks, Spline::QUADRATIC);
		for (auto y : vy) fouty << y << "\t";
		fouty << endl;

	}

	S.print();

	fout1.open("iebt2d_u1.txt");
	S.species_vec[0]->printCohortVector(fout1);
	fout1.close();

	fout.close();
	foutx.close();
	fouty.close();

	// cout << S.u0_out(S.current_time)[0] << endl; 
	// if (abs(S.u0_out(S.current_time)[0] - 1.30639) < 1e-5) return 0;
	// else return 1;

}

