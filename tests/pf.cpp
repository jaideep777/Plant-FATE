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

int main(int argc, char ** argv){

	string pfile = "tests/params/p.ini";
	if (argc > 1) pfile = argv[1];

	for (int i=0; i<1; ++i){
		Simulator sim(pfile);
		sim.init(1000, 1010);
		sim.simulate();
		sim.close();
	}


	// // 1. eCO2 run for 2 different zetas
	// vector<double> zeta = {0.2, 0.08};

	// for (int i=0; i<1; ++i){
	// 	Simulator sim(pfile);
	// 	sim.traits0.zeta = zeta[i];
	// 	sim.expt_dir = "HIST_ELE_zeta_" + to_string(zeta[i]);
	// 	sim.init(-1000, 5000);
	// 	sim.simulate();
	// 	sim.close();
	// }


	// // Test effect of zeta
	// vector<double> zeta = myseq(0.03, 0.3, 13);

	// std::for_each(
	// 	std::execution::par_unseq,
	// 	zeta.begin(),
	// 	zeta.end(),
	// 	[](double zz){
	// 		Simulator sim("tests/params/p.ini");
	// 		sim.traits0.zeta = zz;
	// 		sim.expt_dir = "HIST_zeta_" + to_string(zz);
	// 		sim.init(1000, 1500);
	// 		sim.simulate();
	// 		sim.close();
	// 	});




	// // Effect of CO2
	// vector<double> co2_vec = myseq(360, 600, 12);
	// for (auto cc : co2_vec){
	// 	Simulator sim("tests/params/p.ini");
	// 	sim.expt_dir = "scan_co2_" + to_string(cc);
	// 	sim.E.clim.co2 = cc;
	// 	sim.init(1000, 1500);
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


