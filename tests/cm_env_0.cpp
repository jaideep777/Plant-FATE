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

	Solver S(SOLVER_CM);
	S.setEnvironment(&E);
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 4);
	S.print();
	
	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;

	if (fabs(E.evalEnv(0,0) - 0.3821924) > 1e-6) return 1;
	
	return 0;
	
}

