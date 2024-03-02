#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms_precompute.h"

int main(){

	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_CM);
	S.control.cm_use_log_densities = true;
	S.control.cm_grad_dx = 0.001;
	S.control.max_cohorts = 26;
	S.setEnvironment(&E);
	S.addSpecies(25, 0, 1, false, &spp, 4, 2);
	S.resetState();
	S.initialize();
	S.species_vec[0]->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	ofstream fout("cm_testmodel.txt");

	fout << S.current_time << "\t" << 0 << "\t";
	for (auto y : S.state){ fout << y << "\t";} fout << "\n";
	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		//cout << S.current_time << "\t" << S.species_vec[0]->xsize() << " " << S.u0_out()[0] << "\t" << S.species_vec[0]->get_boundary_u() << "\n";
		//cout << S.u0_out() << "\n";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << "Number of calls to p/g/m/f (static) : " << Cohort<TestModel>::np << " " << Cohort<TestModel>::ng << " " << Cohort<TestModel>::nm << " " << Cohort<TestModel>::nf << endl;
	cout << "Number of calls to derivs           : " << S.odeStepper.get_fn_evals() << endl;

	cout << S.u0_out(S.current_time)[0] << endl;
	// test value is from R code	
	//if (abs(S.u0_out()[0] - 1.556967) < 1e-5) return 0;  // this is when integrate_x BC is not included
	if (abs(S.u0_out(S.current_time)[0] - 1.397015) < 1e-5) return 0;  // this is when integrate_x BC IS included

	else return 1;
  
}

