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

	Solver S(SOLVER_EBT);
	S.control.ebt_grad_dx = 0.001;
	//S.control.ode_method = "rk4";
	//S.control.ode_rk4_stepsize = 0.01;
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

	vector<double> breaks = myseq(0,1,26);

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";

		//cout << "HERE" << endl;	
		vector<double> v = S.getDensitySpecies(0, breaks);
				
		cout << S.current_time << " " << S.species_vec[0]->xsize() << " " << S.u0_out(t)[0] << endl;
		//for (auto y : S.state) fout << y << "\t";
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	S.print();	
	cout << S.u0_out(S.current_time)[0] << endl;
	if (abs(S.u0_out(S.current_time)[0]-1.436407) < 2e-5) return 0;
	else return 1;
}

