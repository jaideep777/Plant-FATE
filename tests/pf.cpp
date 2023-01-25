#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include "plantfate.h"

using namespace std;

int main(){
	
	{
		Simulator sim("tests/params/p_ele_base.ini");
		sim.init(1000, 3000);
		sim.simulate();
		sim.close();
	}

	{
		Simulator sim("tests/params/p_ele_hi.ini");
		sim.init(1000, 3000);
		sim.simulate();
		sim.close();
	}

	{
		Simulator sim("tests/params/p_ele_low.ini");
		sim.init(1000, 3000);
		sim.simulate();
		sim.close();
	}


}


