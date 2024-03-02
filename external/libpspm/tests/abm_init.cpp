#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	//TestModel M;
	Environment E;
	Species<TestModel> spp;

	Solver S(SOLVER_ABM);
	S.control.abm_n0 = 1000;
	S.setEnvironment(&E);
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 4);

	ofstream fout("abm_init.txt");
	S.species_vec[0]->printCohortVector(fout);
	fout.close();

	// S.print();
	
	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;

	if (fabs(E.evalEnv(0,0) - 0.3821924) > 5e-3) return 1;
	
	return 0;
	
}

