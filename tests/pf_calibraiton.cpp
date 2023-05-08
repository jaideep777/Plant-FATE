#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include "plantfate.h"

using namespace std;

int main(){

	Simulator sim("tests/params/p_calib.ini");
	//sim.expt_dir = "AMB_HD_CALIB";
	sim.S.control.ebt_merge_dxcut = 0.001;
	sim.E.clim.co2 = 414;
	sim.init(1000, 2000);
	sim.simulate();
	sim.close();

}


