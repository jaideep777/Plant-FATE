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

	Solver S(SOLVER_IFMU);
	S.setEnvironment(&E);
	// S.control.ode_ifmu_stepsize = 0.001;
	S.control.ifmu_centered_grids = true; // FIXME: 'false' seems more correct, but to cross check with 1D impl, set to true. Investigate
	S.addSpecies({30}, {0}, {1}, {false}, &spp, 4, -1);
	// S.species_vec[0]->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	// S.resetState();
	// S.initialize();
	S.print();
	
	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;
	
	ofstream fout("ifmu_testmodel_equil.txt");

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << " " << S.u0_out(t)[0] << "\n";
		
		vector<double> breaks = myseq(0,1,26);
		vector<double> v = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << setprecision(7) << S.u0_out(S.current_time)[0] << endl; 
	// if (abs(S.u0_out(S.current_time)[0] - 1.306386) < 5e-4) return 0; // this value is expected if xb = X[0] (probably not correct)
	if (abs(S.u0_out(S.current_time)[0] - 1.306023) < 1e-6) return 0; // this value is expected if xb = x[0] (probably correct)
	else return 1;

}

