#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>
#include <execution>

#include "plantfate.h"

using namespace std;

int main(){

	// Test effect of zeta
	vector<double> zeta = myseq(0.03, 0.3, 13);

	std::for_each(
		std::execution::par_unseq,
		zeta.begin(),
		zeta.end(),
		[](double zz){
			Simulator sim("tests/params/p.ini");
			sim.traits0.zeta = zz;
			sim.expt_dir = "HIST_zeta_" + to_string(zz);
			sim.init(1000, 1500);
			sim.simulate();
			sim.close();
		});



	// // Effect of zeta x CO2
	// vector<double> zeta = {0.08, 0.2};

	// for (int i=0; i<2; ++i){
	// 	Simulator sim("tests/params/p.ini");
	// 	sim.traits0.zeta = zeta[i];
	// 	sim.expt_dir = "HIST_ELE_zeta_" + to_string(zeta[i]);
	// 	sim.init(1000, 3000);
	// 	sim.simulate();
	// 	sim.close();
	// }



	// Ensure that 1 mo timestep is good enough... 
	// {
	// 	Simulator sim("tests/params/p.ini");
	// 	sim.timestep = 1.0/12.0/2.0;
	// 	sim.expt_dir = "BASE_dt_2wk";
	// 	sim.init(1000, 2000);
	// 	sim.simulate();
	// 	sim.close();
	// }

	// {
	// 	Simulator sim("tests/params/p.ini");
	// 	sim.timestep = 1.0/12.0;
	// 	sim.expt_dir = "BASE_dt_4wk";
	// 	sim.init(1000, 2000);
	// 	sim.simulate();
	// 	sim.close();
	// }



	// Test effect of eCO2 under different zeta
	// {
	// 	Simulator sim("tests/params/p_ele_base.ini");
	// 	sim.init(1000, 2000);
	// 	sim.simulate();
	// 	sim.close();
	// }

	// {
	// 	Simulator sim("tests/params/p_ele_hi.ini");
	// 	sim.init(1000, 3000);
	// 	sim.simulate();
	// 	sim.close();
	// }

	// {
	// 	Simulator sim("tests/params/p_ele_low.ini");
	// 	sim.init(1000, 3000);
	// 	sim.simulate();
	// 	sim.close();
	// }


}


