#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include "plantfate.h"

using namespace std;

int main(){

	// reference run
	{
		Simulator sim("tests/params/p.ini");
		sim.continuePrevious = false;
		sim.out_dir = sim.parent_dir + "/" + "cont_test_ref";
		sim.init(1000, 1350);
		sim.simulate();
		sim.close();
	}

	// spinup run
	{
		Simulator sim("tests/params/p.ini");
		sim.continuePrevious = false;
		sim.out_dir = sim.parent_dir + "/" + "cont_test_spinup";
		sim.init(1000, 1200);
		sim.simulate();
		sim.close();
	}

	// continuation run
	{
		Simulator sim("tests/params/p.ini");
		sim.continuePrevious = true;
		sim.continueFrom_stateFile    = sim.parent_dir + "/" + "cont_test_spinup/pf_saved_state.txt"; 
		sim.continueFrom_configFile   = sim.parent_dir + "/" + "cont_test_spinup/pf_saved_config.ini"; 
		sim.out_dir = sim.parent_dir + "/" + "cont_test_main";
		sim.init(1000, 1350);
		sim.simulate();
		sim.close();
	}


}

