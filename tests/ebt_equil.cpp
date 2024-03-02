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

	// TestModel M;
	// Species<TestModel> spp(M);
	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_EBT);
	S.control.ebt_grad_dx = 0.001;
	S.control.cohort_insertion_dt = 0.05;
	S.setEnvironment(&E);
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 4, -1);

	S.print();
	
	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;
	if (fabs(E.evalEnv(0,0) - 0.3802) > 1e-5) return 1;

	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	
	ofstream fout("ebt_testmodel_equil.txt");

	vector<double> breaks = myseq(0,1,26);

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";

		vector<double> breaks = myseq(0,1,26);
		vector<double> v = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	S.print();	
	cout << setprecision(6) << S.u0_out(S.current_time)[0] << endl;
	if (abs(S.u0_out(S.current_time)[0]-0.97504) < 1e-6) return 0;
	else return 1;
}

