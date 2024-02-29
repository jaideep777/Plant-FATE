#ifndef PLANT_FATE_LIGHT_ENVIRONMENT_H_
#define PLANT_FATE_LIGHT_ENVIRONMENT_H_

#include <iostream>
#include <vector>
#include <environment_base.h>
#include <solver.h>

namespace pfate{
namespace env{

class LightEnvironment : public EnvironmentBase {
	public:
	bool use_ppa = true;

	// Smooth environment 

	// PPA environment	
	int n_layers;
	double total_crown_area;
	std::vector<double> z_star;
	std::vector<double> fapar_tot;
	std::vector<double> canopy_openness;
	
	public:
	LightEnvironment();

	double projected_crown_area_above_z(double t, double z, Solver *S);
	double fapar_layer(double t, int layer, Solver *S);

	void computeEnv(double t, Solver * sol, std::vector<double>::iterator S, std::vector<double>::iterator dSdt) override;

	virtual void print(double t);
	
};

} // env
} // namespace pfate

/*#include "../src/light_environment.tpp"*/

#endif

