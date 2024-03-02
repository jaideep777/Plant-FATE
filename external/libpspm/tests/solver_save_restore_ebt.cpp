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

	double t;
	ofstream fout("ebt_testmodel_equil.txt");
	{
		Solver S(SOLVER_EBT);
		S.control.ebt_grad_dx = 0.001;
		S.control.sync_cohort_insertion = true;
		//S.control.ode_method = "rk4";
		//S.control.ode_rk4_stepsize = 0.01;
		S.setEnvironment(&E);
		S.addSpecies({25}, {0}, {1}, {false}, &spp, 4, -1);
		S.print();
		S.print();
		
		E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
		cout << E.evalEnv(0,0) << endl;

		S.print();
		//for (auto s : S.state) cout << s << " "; cout << endl;
		
		vector<double> breaks = myseq(0,1,26);
	
		for (t=0.05; t <= 4; t=t+0.05) {
			S.step_to(t);
			fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";

			vector<double> breaks = myseq(0,1,26);
			vector<double> v = S.getDensitySpecies1D(0,0, breaks, Spline::QUADRATIC);
			for (auto y : v) fout << y << "\t";
			fout << endl;
		}

		ofstream fouts("solver_state_save.txt");
		S.save(fouts);
		fouts.close();

		cout << "  ------------ compare saved and restored solver -------------\n";
		S.print();
	}

	// What follows is like a new instance of a solver that's initialized from file

	{
		ifstream fins("solver_state_save.txt");
		Solver S(SOLVER_EBT);
		S.control.ebt_grad_dx = 0.001;
		S.control.sync_cohort_insertion = true;
		S.setEnvironment(&E);

		Species<TestModel> spp_proto;
		S.restore(fins, {&spp_proto});
		S.print();

		for (; t <= 8; t=t+0.05) {
			S.step_to(t);
			fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";

			vector<double> breaks = myseq(0,1,26);
			vector<double> v = S.getDensitySpecies1D(0,0, breaks, Spline::QUADRATIC);
			for (auto y : v) fout << y << "\t";
			fout << endl;
		}

		fout.close();

		S.print();	
		cout << setprecision(10) << S.u0_out(S.current_time)[0] << endl;
		if (abs(S.u0_out(S.current_time)[0]-0.9750400229) < 1e-7) return 0;
		else return 1;
	}

}

