#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_FMU, "lsoda");
	S.control.ode_eps = 1e-4;
	S.setEnvironment(&E);
	S.addSpecies(25, 0, 1, false, &spp, 4, 2);
	S.species_vec[0]->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.resetState();
	S.initialize();
	S.print();
	
	ofstream fout("fmu_testmodel.txt");

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << " " << S.u0_out(t)[0] << "\n";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << S.u0_out(S.current_time)[0] << endl; 
	cout << "Number of fn evaluations = " << S.odeStepper.get_fn_evals() << "\n";
	if (abs(S.u0_out(S.current_time)[0] - 1.468232) < 1e-5) return 0;
	else return 1;

}

