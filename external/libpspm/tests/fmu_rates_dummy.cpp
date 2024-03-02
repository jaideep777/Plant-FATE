#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_FMU);
	S.setEnvironment(&E);
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 4, -1);
	S.species_vec[0]->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.print();
	
	E.computeEnv(0,&S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;

	S.calcRates_FMU(1, S.state.begin(), S.rates.begin());  // dummy rates calc rates(X=X0, U=U0, t=1, E=E(U0))
	//S.print();
	//	S.step_to(1);

	// This is with input birth flux set to -1
	vector <double> rates_exp = {
		-0.355862928, -0.450926505, -0.417535387, -0.354046196, 
		-0.301067280, -0.256580544, -0.219012561, -0.187124217,
	    -0.159930669, -0.136642529, -0.116622139, -0.099350697,
	    -0.084403290, -0.071429752, -0.060139860, -0.050291781,
	    -0.041683001, -0.034143136, -0.027528209, -0.021716055,
	    -0.016602617, -0.012098936, -0.008128703, -0.004910899,
	    -0.001250216};
	for (int i=0; i< 25-2; ++i){ // FIXME: skip last two because depends on first order or second order approx used
		cout << S.rates[i] << " " << rates_exp[i] << endl;
		if ( fabs(S.rates[i] - rates_exp[i]) > 1e-5) return 1;
	}
	
	return 0;
	
}


