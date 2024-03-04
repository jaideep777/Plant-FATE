#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>

#include "plantfate_patch.h"

using namespace std;

int is_equal(const vector<double>& v1, const vector<double>& v2, double tol=1e-6){
	bool b = true;
	for (int i=0; i<v1.size(); ++i){
		cout << "Comparing v["<<i<<"]: " << v1[i] << " " << v2[i] << '\n';
		b = b & (fabs(v1[i]-v2[i])< tol);
	}
	return b? 0:1;
}

int main(int argc, char ** argv){

	string pfile = "tests/params/p_test_v2_evol2.ini";

	pfate::Patch sim(pfile);
	sim.init(2000, 2020);

	return 0;
}


