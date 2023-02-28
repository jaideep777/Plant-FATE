#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include "plantfate.h"

using namespace std;

int main(){

	Simulator sim("tests/params/p.ini");
	sim.expt_dir = "HIST_HD_CALIB";
	sim.init(1000, 3000);
	sim.simulate();
	sim.close();

}


