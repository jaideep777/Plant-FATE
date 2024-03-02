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
	S.addSpecies({{0,.1,.1+1e-7,.1+2e-7, .2, .2+1e-7}}, &s1, 4, 2);
	S.addSpecies({{0,.1,.2,.3,.3+1e-6,.3+2e-6,.3+3e-6, .4, .4+1e-7,.4+1e-6}}, &s2, 4, 2);
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;
	
	auto spp1 = S.species_vec[0];	
	auto spp2 = S.species_vec[1];	


	S.print();

	S.mergeCohorts_EBT();

	S.print();

	return 0;
	
}

