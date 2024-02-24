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

	string pfile = "tests/params/p_base.ini";
	if (argc > 1) pfile = argv[1];

	for (int i=0; i<1; ++i){
		Simulator sim(pfile);
		sim.init(-1000, 5000);

		sim.set_co2File("");

		sim.S.print();
		cout << setprecision(6) << "t\ttc\tCO2414\tGPP414\tCO2614\tGPP614\tpc\n";
		for (double t = 2000; t < 2016; t += 0.5/12.0){
			sim.S.current_time = t;

			sim.S.updateEnv(sim.S.current_time,  sim.S.state.begin(), sim.S.rates.begin());
			sim.E.clim_inst.co2 = 414;
			double co21 = sim.E.clim_inst.co2;
			sim.S.calcRates_EBT(sim.S.current_time, sim.S.state.begin(), sim.S.rates.begin());
			sim.props.update(sim.S.current_time, sim.S);
			double gpp1 = sim.props.gpp*0.5;

			sim.S.updateEnv(sim.S.current_time,  sim.S.state.begin(), sim.S.rates.begin());
			sim.E.clim_inst.co2 = 614;
			double co22 = sim.E.clim_inst.co2;
			sim.S.calcRates_EBT(sim.S.current_time, sim.S.state.begin(), sim.S.rates.begin());
			sim.props.update(sim.S.current_time, sim.S);
			double gpp2 = sim.props.gpp*0.5;

			cout << t << "\t" << sim.E.clim_inst.tc << "\t" << co21 << "\t" << gpp1 << "\t" << co22 << "\t" << gpp2 << "\t" << (gpp2/gpp1-1)*100 << "\n";

		}


	}


}


