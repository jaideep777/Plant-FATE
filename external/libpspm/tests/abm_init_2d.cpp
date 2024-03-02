#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "wave_model_2d.h"

int main(){

	//TestModel M;
	WaveEnv E;
	Species<Wave> spp;

	Solver S(SOLVER_ABM);
	S.control.abm_n0 = 20000;
	S.setEnvironment(&E);
	S.addSpecies({100, 100}, {0, 0}, {10,10}, {false, false}, &spp, 0, -1);

	ofstream fout("abm_init_2d.txt");
	S.species_vec[0]->printCohortVector(fout);
	fout.close();

	// S.print();
	
	// E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	// cout << E.evalEnv(0,0) << endl;

//	if (fabs(E.evalEnv(0,0) - 0.3821924) > 1e-6) return 1;
	
	return 0;
	
}

