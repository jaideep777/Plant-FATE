#include <iostream>
#include <solver.h>

#include "test_model_2_ms.h"

int main(){
	
	//TestModel M;
	Environment E;
	Species<TestModel> spp1;
	Species<TestModel> spp2;
	Species<TestModel> spp3;

	Solver S(SOLVER_CM);
	S.setEnvironment(&E);
	S.addSpecies({10}, {0}, {1}, {false}, &spp1, 4, 10);
	S.addSpecies({5}, {0}, {0.5}, {false}, &spp2, 4, 5);
	S.addSpecies({2}, {0}, {1}, {false}, &spp3, 4, 2);
	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	S.print();

	S.addCohort_CM();
	S.print();

	S.addCohort_CM();
	S.print();

	// test user inserting a user-defined cohort
	Cohort<TestModel> bc;
	bc.x = {1.555};
	bc.u = 900;
	bc.mortality = 1;
	bc.viable_seeds = 2;
	bc.heart_mass = 200;
	bc.sap_mass = 300;
	spp1.addCohort(bc);
	S.print();
	
	//S.removeCohort_CM();
	//S.print();

	Species<TestModel> spp4;
	Solver S1(SOLVER_CM);
	S1.setEnvironment(&E);
	S1.control.max_cohorts = 5;
	S1.addSpecies({{.1,.4,.5,.6,1.0}}, &spp4, 4, 10);
	S1.print();

	S1.removeCohort_CM();
	S1.print();
		
	return 0;
}

