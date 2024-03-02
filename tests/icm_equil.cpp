#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

int main(){

	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_ICM);
	S.setEnvironment(&E);
	S.control.cm_use_log_densities = false;
	S.control.cm_grad_dx = {0.001};
	S.control.ode_ifmu_stepsize = 0.001;
	S.control.max_cohorts = 26;
	S.control.update_cohorts = true;
	S.control.cm_remove_cohorts = true;
	S.control.cohort_insertion_dt = 0.05;
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 4, -1);
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;

	ofstream fout("icm_testmodel_equil.txt");

	//fout << S.current_time << "\t" << 0 << "\t";
	//for (auto y : S.state){fout << y << "\t";} fout << "\n";
	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << "\t" << S.species_vec[0]->xsize() << "\t" << S.u0_out(t)[0] << "\n";
		//cout << S.current_time << "\t" << S.species_vec[0]->xsize() << " " << S.u0_out()[0] << "\t" << S.species_vec[0]->get_boundary_u() << "\n";
		//cout << S.u0_out() << "\n";

		vector<double> breaks = myseq(0,1,26);
		vector<double> v = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	// S.print();

	cout << setprecision(10) << S.u0_out(S.current_time)[0] << endl;
	// test value is from R code	
	//if (abs(S.u0_out()[0] - 1.556967) < 1e-5) return 0;  // this is when integrate_x BC is not included
	if (abs(S.u0_out(S.current_time)[0] - 0.998272) < 1e-6) return 0;  // this is when integrate_x BC IS included
	else return 1;
  
}

