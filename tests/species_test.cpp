#include <iostream>
#include <vector>
#include <cohort.h>
#include <species.h>
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
int check(const vector<T>& v1, const vector<T>& v2){
	if (v1.size() != v2.size()) return 1;

	for (int i=0; i<v1.size(); ++i){
		if (fabs(v1[i]-v2[i]) > 1e-6) return 1;
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

	cout << "Sort cohorts by body size\n-----------------------\n";
	Si.sortCohortsDescending(0);
	Si.print();
	nerrors += check_species_states(Si, {20, 0.01, 20, 10, 0.5, 5});

	cout << "Sort cohorts by energy reserve\n-----------------------\n";
	Si.sortCohortsDescending(1);
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 20, 0.01});

	// Now insert a EBT boundary cohort and check that it can be skipped in sorting
	Cohort<Insect> Cb; Cb.set_size({100, 0.02}); Cb.u = 0.1;
	Si.addCohort(Cb);
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 20, 0.01, 100, 0.02});

	cout << "Sort cohorts by body size, exclude BC\n-----------------------\n";
	Si.sortCohortsDescending(0, 1);
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 20, 0.01, 0.5, 5, 100, 0.02});


	cout << "Set state of a specific cohort\n-----------------------\n";
	Si.setX(1, {50, 20});
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 50, 20, 0.5, 5, 100, 0.02});
	
	auto g1 = Si.growthRate(2, 1, &E); // cohort 2 {0.5, 5}
	cout << g1 << '\n';
	{
	nerrors += check(g1, {0.5*5*0.1, 1*0.5*0.1});
	nerrors += check(Si.mortalityRate(2, 1, &E), 0.5/5);
	nerrors += check(Si.birthRate(2, 1, &E), 0.05);
	}

	cout << "Testing removal of a specific cohort\n-----------------------\n";
	Si.markCohortForRemoval(1);
	Si.removeMarkedCohorts();
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 100, 0.02});

	cout << "Testing removal of a dead cohorts\n-----------------------\n";
	Cohort<Insect> Cdead(I1); Cdead.set_size({1000, 1e-4}); Cdead.u = 1e-15;
	Si.addCohort(Cdead);
	Si.addCohort(Cdead);
	Si.addCohort(Cdead);
	Si.print();

	Si.markDeadCohorts(1e-10); // Note: This will not remove cohort at index J-1. We must ensure during cohort sorting that J-1 is indeed the boundary cohort
	Si.removeMarkedCohorts(); // Note: This will not remove cohort at index J-1. We must ensure during cohort sorting that J-1 is indeed the boundary cohort
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 100, 0.02, 1000, 1e-4});


	cout << "Testing addCohort()\n-----------------------\n";
	Si.addCohort();
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 100, 0.02, 1000, 1e-4, 0.5, 0.01});


	// nerrors += check(Sp.growthRate(2, 1, &E), {32});


	cout << "Testing mortalityRateGradient()\n-----------------------\n";
	vector<double> m_mx = Si.mortalityRateGradient(0, 1, &E, {1e-6, 1e-6});
	cout << "Mort rate gradient = " << m_mx << '\n';
	{
	auto& c = Si.getCohort(0);
	nerrors += check(m_mx, {c.w/c.e, 1/c.e, -c.w/c.e/c.e});
	//                        ^ This is analytical gradient calculation
	}


	// JJ FIXME: Double check this calculation.
	cout << "Testing growthRateGradient()\n-----------------------\n";
	vector<vector<double>> g_gx = Si.growthRateGradient(0, 1, &E, {1e-6, 1e-5});
	cout << "Growth rate gradient:\n";
	for (int i=0; i<Si.istate_size+1; ++i){
		cout << vector<std::string>{"   g:"," gx1:", " gx2:"}[i] << " ";
		for (int j=0; j<Si.istate_size; ++j){
			cout << setw(10) << g_gx[i][j] << ' ';
		}
		if (i==0) cout << "\n----------------------------";
		cout << '\n';
	}
	nerrors += check(g_gx[0], {20, 2});
	{
	auto& c = Si.getCohort(0);
	nerrors += check(g_gx[1], {c.e*0.1, E.E*0.1});
	nerrors += check(g_gx[2], {c.w*0.1, 0});
	//                        ^ This is analytical gradient calculation
	}


	cout << "Testing accumulators\n-----------------------\n";
	Sp.print();
	Sp.initAccumulators(1, &E);
	Sp.print();
	vector<double> acc(4, 0);
	auto it = acc.begin();
	Sp.copyAccumulatorsToState(it);
	nerrors += check(acc, vector<double>(4, 2.03));

	acc = {1.01, 2.02, 3.03, 4.04};
	it = acc.begin();
	Sp.copyAccumulatorsToCohorts(it);
	Sp.print();
	nerrors += check(Sp.getCohort(0).root_mass, 1.01);
	nerrors += check(Sp.getCohort(1).root_mass, 2.02);
	nerrors += check(Sp.getCohort(2).root_mass, 3.03);
	nerrors += check(Sp.getCohort(3).root_mass, 4.04);

	it = acc.begin();
	Sp.accumulatorRates(it);
	cout << "Accum var rates: " << acc << '\n';
	nerrors += check(acc, {-0.001, -0.102, -0.203, -0.304});


	cout << "Testing density initialization\n-----------------------\n";
	for (int i=0; i<Si.xsize(); ++i) Si.setU(i, Si.init_density(i, &E));
	Si.print();	
	vector<double> expected_u = {20, 0.25, 0.2, 0.01, 0.0005};
	for (int i=0; i<Si.xsize(); ++i) nerrors += check(Si.getU(i), expected_u[i]);


	cout << "Testing cohort grouping\n-----------------------\n";
	Cohort<Insect> Ci1(I1);
	Ci1.set_size({2,5.5});
	Si.addCohort(Ci1);
	Si.addCohort(Ci1);
	Si.addCohort(Ci1);
	Ci1.set_size({0.5, 0.01});
	Si.addCohort(Ci1);

	for (int i=0; i<Si.xsize(); ++i) Si.set_birthTime(i, (Si.xsize()-i)*0.1);

	Si.save(cout);
	Si.markDenseCohorts(1e-6);
	Si.save(cout);
	Si.removeMarkedCohorts();
	Si.save(cout);

	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 100, 0.02, 1000, 0.0001, 2, 5.5, 2, 5.5, 0.5, 0.01});


	cout << "Testing save / restore\n-----------------------\n";
	ofstream fout("Si_save.txt");
	Si.save(fout);
	fout.close();

	ifstream fin("Si_save.txt");
	Species<Insect> Si_restored(I1);
	Si_restored.restore(fin);
	Si_restored.save(cout);
	// Verify manually that console output is same as in the saved file.


	if (nerrors == 0) cout << "******* ALL TESTS PASSED ***********\n";
	else cout << "xxxxxx " << nerrors << " TESTS FAILED xxxxxxxxxxx\n";
	
	return nerrors;

// 	Species<Plant> Sv(array<double,1> {1,2,3,4,5});
// 	Sv.setX(0, 1.5);
// 	Sv.setX(2, 3.5);
// 	Sv.print();

// 	Species<Plant> S(P);  // create species S with prototype P

// 	Species<Insect> I(array<double,1> {1.1,2.1,3.1});

// 	Species_Base * S1 = &S;
// 	//S1->print();
// 	cout << "S size: " << S1->get_maxSize() << "\n";

// 	Species_Base * S2 = &I;
// 	//S2->print();
// 	cout << "I size: " << S2->get_maxSize() << "\n";

// 	LightEnv env;

// 	Solver sol(SOLVER_EBT);
// 	sol.setEnvironment(&env);
// 	sol.addSpecies(array<double,1> {1,2,3,4,5}, &S, 1, 1);
// 	sol.addSpecies(array<double,1> {1.1,2.1,3.1}, &I, 0, 1);
// 	sol.resetState();
// 	sol.initialize();
// 	//sol.addCohort_EBT();
// 	sol.print();	

// 	I.setX(2, 0.05);
// 	I.setU(2, 0.2);
// 	sol.print();	
// 	cout << "Insect BF = " << sol.calcSpeciesBirthFlux(1,0) << "\n";
	
// //	sol.preComputeSpecies(1,0);
// 	sol.print();	
// 	cout << "Insect BF (after precompute) = " << sol.calcSpeciesBirthFlux(1,0) << "\n";

// 	Cohort<Plant> C;
// 	C.print(); cout << "\n";
// 	C.set_size(10);
// 	C.print(); cout << "\n";
// 	cout << C.growthRate(C.height, 0, sol.env) << "\n";

}



