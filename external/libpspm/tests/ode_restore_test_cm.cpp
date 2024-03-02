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

	Solver S(SOLVER_CM);
	S.setEnvironment(&E);
	S.control.cm_use_log_densities = false;
	S.control.cm_grad_dx = {0.001};
	S.control.max_cohorts = 26;
	S.control.cm_remove_cohorts = true;
	S.control.sync_cohort_insertion = true;
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 4, -1);
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;

	ofstream fout("cm_testmodel_ode_restored.txt");

	// fout << S.current_time << "\t" << 0 << "\t";
	// for (auto y : S.state){ fout << y << "\t";} fout << "\n";
	for (double t=0.05; t <= 4; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << "\t" << S.species_vec[0]->xsize() << " " << S.u0_out(S.current_time)[0] << "\t" << S.species_vec[0]->get_boundary_u() << "\n";
		//cout << S.u0_out() << "\n";

		vector<double> breaks = myseq(0,1,26);
		vector<double> v = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	// S.print();
	S.odeStepper.save(cout); 
	cout << '\n';

	ofstream fouts("ode_state_save.txt");
	S.odeStepper.save(fouts);
	fouts.close();
	S.odeStepper.reset(0,1e-6,1e-6);
	ifstream fins("ode_state_save.txt");
	S.odeStepper.restore(fins);

	S.odeStepper.save(cout);
	cout << '\n';

	for (double t=4.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << "\t" << S.species_vec[0]->xsize() << " " << S.u0_out(S.current_time)[0] << "\t" << S.species_vec[0]->get_boundary_u() << "\n";
		//cout << S.u0_out() << "\n";

		vector<double> breaks = myseq(0,1,26);
		vector<double> v = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << "Number of calls to p/g/m/f (static) : " << Cohort<TestModel>::np << " " << Cohort<TestModel>::ng << " " << Cohort<TestModel>::nm << " " << Cohort<TestModel>::nf << endl;
	cout << "Number of calls to derivs           : " << S.odeStepper.get_fn_evals() << endl;

	cout << setprecision(10) << S.u0_out(S.current_time)[0] << endl;
	// test value is from R code	
	//if (abs(S.u0_out()[0] - 1.556967) < 1e-5) return 0;  // this is when integrate_x BC is not included
	if (abs(S.u0_out(S.current_time)[0] - 0.976177) < 1e-5) return 0;  // this is when integrate_x BC IS included

	else return 1;
  
}

