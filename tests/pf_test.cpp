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

int is_equal(const vector<double>& v1, const vector<double>& v2, double tol=1e-6){
	bool b = true;
	for (int i=0; i<v1.size(); ++i){
		cout << "Comparing v["<<i<<"]: " << v1[i] << " " << v2[i] << '\n';
		b = b & (fabs(v1[i]-v2[i])< tol);
	}
	return b? 0:1;
}

int main(int argc, char ** argv){

	string pfile = "tests/params/p_test.ini";
	if (argc > 1) pfile = argv[1];

	int err = 0;

	for (int i=0; i<1; ++i){
		Simulator sim(pfile);
		// sim.expt_dir = sim.expt_dir + "_414ppm";
		sim.E.clim.co2 = 414;
		sim.init(0, 100);
		sim.simulate();

		vector<double> ba = sim.cwm.ba_vec;
		for (auto& b : ba) b*=1e4;
		cout << setprecision(10) << "Basal areas [m2/Ha]: " << ba << '\n';

		// err = is_equal(ba, {1.568911254, 7.061621143, 27.19173438});  // this was the output at yearly step_to when using 1d libpspm, reproduced in 2d
		// err = is_equal(ba, {1.569239985, 7.069009343, 27.21039394});  // this is the output  at yearly step_to after ODE bugfix in libpspm @8d0eb68
		// err = is_equal(ba, {1.569239985, 7.069009343, 27.21039394});  // this is the output  at yearly step_to after ODE bugfix in libpspm @118d623
		// err = is_equal(ba, {1.569460052, 7.060768412, 27.17424297});     // this is the output  at BIWEEKLY step_to. This difference is due to insertion of an extra point in R0 moving averager at integer t in yearly stepping - this doesnt happen at biweekly stepping, which seems more correct. This is the best level of investigattion I can do for now, so I'm going to take this as the expected test result and proceed with other stuff.
		err = is_equal(ba, {1.569063555, 7.059273716, 27.16867421});     // this is the output  at BIWEEKLY step_to after moving climate update to before step. Difference is because not climate is not updated in after_step() before computation of seed rain

		sim.close();
	}

	return err;

	// for (int i=0; i<1; ++i){
	// 	Simulator sim(pfile);
	// 	sim.expt_dir = sim.expt_dir + "_614ppm";
	// 	sim.E.clim.co2 = 614;
	// 	sim.init(1000, 5000);
	// 	sim.simulate();
	// 	sim.close();
	// }


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

