#ifndef PLANT_FATE_PSPM_INTERFACE_H_
#define PLANT_FATE_PSPM_INTERFACE_H_

#include <solver.h>
#include "light_environment.h"
#include "climate.h"
#include "plant.h"

class PSPM_Dynamic_Environment : public EnvironmentBase, public env::LightEnvironment, public env::Climate{
	public:
	
	void computeEnv(double t, Solver *S){
		LightEnvironment::computeEnv(t, S);
	}
	
	void print(double t){
		Climate::print(t);
		LightEnvironment::print();
	}
};


class PSPM_Plant : public plant::Plant {
	public:
	
	double t_birth = 0;

	std::vector<std::string> varnames = {"|lma|", "D", "g", "mort", "lai", "seeds"}; // header corresponding to the print function below
	std::vector<std::string> statevarnames = {"lai", "mort", "seedpool"};                 // header corresponding to state output

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	PSPM_Plant(); 

	void set_size(double _x);

	double init_density(double x, void * _env, double input_seed_rain);

	void preCompute(double x, double t, void * _env);

	double establishmentProbability(double t, void * _env);

	double growthRate(double x, double t, void * _env);
	double mortalityRate(double x, double t, void * _env);
	double birthRate(double x, double t, void * _env);
	
	void init_state(double t, void * _env);

	std::vector<double>::iterator set_state(std::vector<double>::iterator &it);

	std::vector<double>::iterator get_state(std::vector<double>::iterator &it);
	std::vector<double>::iterator get_rates(std::vector<double>::iterator &it);

	void print(std::ostream &out = std::cout);

};



#endif
