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

	sim.init();

	sim.simulate();

	sim.close();
}

