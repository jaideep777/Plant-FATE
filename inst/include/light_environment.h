#ifndef PLANT_FATE_ENV_LIGHT_ENVIRONMENT_H_
#define PLANT_FATE_ENV_LIGHT_ENVIRONMENT_H_

#include <iostream>
#include <vector>

namespace env{

class LightEnvironment {
	public:
	bool use_ppa = false;

	// Smooth environment 

	// PPA environment	
	int n_layers;
	double total_crown_area;
	std::vector<double> z_star;
	std::vector<double> fapar_tot;
	std::vector<double> canopy_openness;
	
	
	public:
	LightEnvironment();
	void print();
	
};

} // env

/*#include "../src/light_environment.tpp"*/

#endif

