#include <iostream>
#include <cmath>
#include <vector>

#include "light_environment.h"

namespace env{

LightEnvironment::LightEnvironment(){ 
	n_layers = 0;
	z_star = {0};
	canopy_openness = {1};
	z_star.reserve(20);
	canopy_openness.reserve(20);
}


void LightEnvironment::print(){
	std::cout << "PPA:" << "\n";
	std::cout << "z* (" << n_layers << ") = "; for (auto z:z_star) std::cout << z << " "; 
	std::cout << "\ncanopy openness* = "; for (auto z:canopy_openness) std::cout << z << " "; 
	std::cout << "\n";
	//std::cout << "Light profile Interpolator:" << "\n";
	//light_profile.print();
}


} // env


