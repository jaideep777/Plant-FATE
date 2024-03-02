#include <iostream>
#include <vector>
#include <cohort.h>
using namespace std;

#include "test_model_plant_insect.h"

template <class Model>
class TestSpecies{
	public:
	Cohort<Model> c;
	auto getX(int i){
		return c.x;
	}
};


int main(){

	LightEnv E;

	Plant P;
	P.set_size({25});
	cout << "P: "; P.print(); cout << "\n";

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
	C2.preCompute(1, &E);
	cout << "Insect g/m/f: " << C2.g << " / " << C2.m << " / " << C2.f << '\n';
	C2.print(); cout << '\n';

	C2.save(cout, 0);

	TestSpecies<Insect> Si;
	Si.c = C2;
	cout << "Cohort inside Insect species has state: " << Si.getX(0) << '\n';
	cout << "Cohort inside Insect species has state type: " << type_name<decltype(Si.getX(0))>() << '\n';

	TestSpecies<Plant> Sp;
	Sp.c = C1;
	cout << "Cohort inside Plant species has state: " << Sp.getX(0) << '\n';
	cout << "Cohort inside Plant species has state type: " << type_name<decltype(Sp.getX(0))>() << '\n';



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



