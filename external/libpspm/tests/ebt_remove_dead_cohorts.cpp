#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

template<class T>
bool almostEqual(const vector<T> &a, const vector <T> &b){
	if (a.size() != b.size()) return false;
	return std::equal(a.begin(), a.end(), b.begin(), [](const T& x, const T& y){return abs(x-y)<1e-5;});
}


int main(){
	Species<TestModel> s1;
	Species<TestModel> s2;
	Environment E;

	Solver S(SOLVER_EBT);
	S.setEnvironment(&E);
	S.addSpecies({5}, {0}, {1}, {false}, &s1, 4, 2);
	S.addSpecies({8}, {0}, {1}, {false}, &s2, 4, 2);
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;
	
	auto spp1 = S.species_vec[0];	
	auto spp2 = S.species_vec[1];	

	spp1->setU(1,1e-11); 
	spp1->setU(2,1e-11); 
	spp2->setU(4,1e-11); 
	spp2->setU(5,1e-11); 

	S.print();

	S.removeDeadCohorts_EBT();

	S.print();

	return 0;
	
}

