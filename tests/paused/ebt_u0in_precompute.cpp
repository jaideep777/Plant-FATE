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

	Solver S(SOLVER_EBT);
	S.control.ebt_grad_dx = 0.001;
	S.setEnvironment(&E);
	S.addSpecies(25, 0, 1, false, &spp, 4, 2);
	S.species_vec[0]->set_bfin_is_u0in(true);
	S.resetState();
	S.initialize();
	S.print();
	
	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;

	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	
	ofstream fout("ebt_testmodel.txt");

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
	
		//vector<double> v = S.cohortsToDensity_EBT(x);
		
		cout << S.current_time << " " << S.species_vec[0]->xsize() << " " << S.u0_out(t)[0] << endl;
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();
	
	cout << "Number of calls to p/g/m/f (static) : " << Cohort<TestModel>::np << " " << Cohort<TestModel>::ng << " " << Cohort<TestModel>::nm << " " << Cohort<TestModel>::nf << endl;
	cout << "Number of calls to derivs           : " << S.odeStepper.get_fn_evals() << endl;
	
	cout << S.u0_out(S.current_time)[0] << endl;
	if (abs(S.u0_out(S.current_time)[0]-1.436407) < 2e-5) return 0;
	else return 1;
}

