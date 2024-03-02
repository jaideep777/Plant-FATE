#include <iostream>
#include <vector>
#include <cohort.h>
#include <species.h>
#include <solver.h>
#include <cmath>
using namespace std;

#include "test_model_plant_insect.h"

template<class Model>
int check_species_states(Species<Model>& spp, const vector<double>& expected){
	for (int i=0, k=0; i<spp.xsize(); ++i){
		auto& c = spp.getCohort(i);
		for (auto& x : c.x){
			if (fabs(x - expected[k]) > 1e-6) return 1;
			++k; 
		}
	}
	return 0;
}

template<class T>
int check(const vector<T>& v1, const vector<T>& v2, double tol=1e-6){
	if (v1.size() != v2.size()) return 1;

	for (int i=0; i<v1.size(); ++i){
		// cout << << v1[i] << " " << v2[i] << '\n';
		if (fabs(v1[i]-v2[i]) > tol) return 1;
	}
	return 0;
}

int check(double v1, double v2){
	if (fabs(v1-v2) > 1e-6) return 1;
	else return 0;
}


int main(){

	int nerrors = 0;

	LightEnv E;

	Plant P;
	P.set_size({25});
	cout << "P: "; P.print(); cout << "\n";
	nerrors += check(P.height, 25);
	nerrors += check(P.crown_area, 625.0);

	P.lma = 70;
	Cohort<Plant> C1;
	C1 = Cohort<Plant> (P);
	C1.set_size({30});
	cout << "C1: "; C1.print(); cout << "\n";

	Insect I1;
	I1.set_size({10, 90});
	cout << "I: "; I1.print(); cout << "\n";

	Cohort<Insect> C2(I1);
	cout << "C2: "; C2.print(); cout << "\n";

	C2.set_size({40,60});
	cout << "C2: "; C2.print(); cout << "\n";

	E.computeEnv(1, nullptr, vector<double>().begin(), vector<double>().begin());
	cout << "Env @ t = 1: " << E.E << '\n';
	nerrors += check(E.E, 1);

	C2.preCompute(1, &E);
	cout << "Insect g/m/f: " << C2.g << " / " << C2.m << " / " << C2.f << '\n';
	C2.print(); cout << '\n';

	C2.save(cout, 0);

	cout << "Testing custom cohort addition\n-----------------------\n";

	Species<Insect> Si(I1); // This constructor does not necessarily set x
	Si.set_xb({0.5, 0.01});
	Si.print();

	Cohort<Insect> C3; C3.set_size({20, 0.01}); C3.u = 1.1;
	Cohort<Insect> C4; C4.set_size({20, 10});   C4.u = 1.3;
	Cohort<Insect> C5; C5.set_size({0.5, 5});   C5.u = 2.5;
	Si.addCohort(C3);
	Si.addCohort(C4);
	Si.addCohort(C5);
	Si.print();

	Species<Plant> Sp(P);
	Sp.set_xb({0.35});
	Sp.n_accumulators = 1;
	Sp.print();	

	Cohort<Plant> Cp3; Cp3.set_size({30}); Cp3.u = 5.1;
	Cohort<Plant> Cp4; Cp4.set_size({31}); Cp4.u = 3.1;
	Cohort<Plant> Cp5; Cp5.set_size({32}); Cp5.u = 1.1;
	Cohort<Plant> Cp6; Cp6.set_size({32.1}); Cp6.u = 0.1;

	Sp.addCohort(Cp3);
	Sp.addCohort(Cp4);
	Sp.addCohort(Cp5);
	Sp.addCohort(Cp6);

	Sp.print();
	nerrors += check_species_states(Sp, {30, 31, 32, 32.1});

	Solver sol_fmu("IFMU", "rk45ck");
	sol_fmu.setEnvironment(&E);
	sol_fmu.addSpecies(vector<int>{5, 3}, vector<double>{1, 0.5}, vector<double>{11,6.5}, {false, false}, (Species_Base*)&Si, 0, -1);
	sol_fmu.addSpecies(vector<int>{10}, vector<double>{0}, vector<double>{20}, {false}, (Species_Base*)&Sp, 1, -1);
	sol_fmu.print();

	Solver sol_ebt("EBT", "rk45ck");
	sol_ebt.setEnvironment(&E);
	sol_ebt.addSpecies(vector<int>{5, 3}, vector<double>{1, 0.5}, vector<double>{11,6.5}, {false, false}, (Species_Base*)&Si, 0, -1);
	sol_ebt.addSpecies(vector<int>{10}, vector<double>{0}, vector<double>{20}, {false}, (Species_Base*)&Sp, 1, -1);
	sol_ebt.addSystemVariables({-50, -60, -70});
	sol_ebt.print();

	cout << "State from FMU solver: " << sol_fmu.state << '\n';
	cout << "State from EBT solver: " << sol_ebt.state << '\n';


	vector<double> expected_state_ebt = 
	    {   -50, -60, -70,
			2,         1.5,         1.2,
			2,         3.5,         2.8,
			2,         5.5,         4.4,
			4,         1.5,         2.4,
			4,         3.5,         5.6,
			4,         5.5,         8.8,
			6,         1.5,         3.6,
			6,         3.5,         8.4,
			6,         5.5,        13.2,
			8,         1.5,         4.8,
			8,         3.5,        11.2,
			8,         5.5,        17.6,
			10,         1.5,           6,
			10,         3.5,          14,
			10,         5.5,          22,
			0,           0,           0,
			1,       6.667,   
			3,       2.857,   
			5,       1.818,   
			7,       1.333,   
			9,       1.053,   
			11,      0.8696,   
			13,      0.7407,   
			15,      0.6452,   
			17,      0.5714,   
			19,      0.5128,   
			0,           0,
			2.03,
			2.03,
			2.03,
			2.03,
			2.03,
			2.03,
			2.03,
			2.03,
			2.03,
			2.03,
			2.03
		};

	nerrors += check(sol_ebt.state, expected_state_ebt, 2e-3);

	sol_ebt.state[0] = -52;
	sol_ebt.state[1] = -62;
	sol_ebt.state[2] = -72;

	sol_ebt.copyStateToCohorts(sol_ebt.state.begin());
	sol_ebt.print();
	nerrors += check(sol_ebt.s_state, {-52, -62, -72});

	expected_state_ebt[0] = -52; 
	expected_state_ebt[0] = -62; 
	expected_state_ebt[0] = -72;

	sol_ebt.copyCohortsToState();
	nerrors += check(sol_ebt.s_state, {-52, -62, -72});


	sol_ebt.save(cout);

	if (nerrors == 0) cout << "******* ALL TESTS PASSED ***********\n";
	else cout << "xxxxxx " << nerrors << " TESTS FAILED xxxxxxxxxxx\n";
	
	return nerrors;


}

