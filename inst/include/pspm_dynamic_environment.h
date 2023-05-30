#ifndef PLANT_FATE_PSPM_DYNAMIC_ENVIRONMENT_H_
#define PLANT_FATE_PSPM_DYNAMIC_ENVIRONMENT_H_

#include <solver.h>
#include "light_environment.h"
#include "climate.h"


/// @ingroup  libpspm_interface
/// @brief    Environment class for interfacing with the PSPM Solver
class PSPM_Dynamic_Environment : public env::LightEnvironment, public env::Climate{
	public:
	PSPM_Dynamic_Environment();
	double projected_crown_area_above_z(double t, double z, Solver *S);
	double fapar_layer(double t, int layer, Solver *S);
	void computeEnv(double t, Solver *S, std::vector<double>::iterator _S, std::vector<double>::iterator _dSdt);
	void print(double t);
};


#endif