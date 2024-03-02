#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	Species<TestModel> M;
	
	Environment E;

	Solver S(SOLVER_FMU);
	S.setEnvironment(&E);
	S.addSpecies({25}, {0}, {1}, {false}, &M, 4);
	S.print();
	
	E.computeEnv(0,&S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;

	if (fabs(E.evalEnv(0,0) - 0.380194) > 1e-6) return 1;
	
	return 0;
	
}

