#ifndef PLANT_PSPM_H_
#define PLANT_PSPM_H_

#include "pspm_environment.h"
#include "plant.h"

#include <solver.h>
#include <individual_base.h>


class PSPM_Plant : public plant::Plant, public IndividualBase<1> {
	public:
	
	double t_birth = 0;
	//double input_seed_rain = 1;	

	double viable_seeds;
	double viable_seeds_dt;

	std::vector<std::string> varnames = {"|lma|", "mort", "fec", "vs", "ha", "hm"};
	std::vector<std::string> statevarnames = {"mort", "vs", "ha", "hm"};

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	PSPM_Plant(); 
	void set_size(const std::array <double, 1>& _x);
	double init_density(void * env, double input_seed_rain);
	void preCompute(double t, void * env);
	double establishmentProbability(double t, void * env);
	std::array<double,1> growthRate(double t, void * env);
	double mortalityRate(double t, void * env);
	double birthRate(double t, void * env);
	
	void init_accumulators(double t, void * _env);
	std::vector<double>::iterator set_accumulators(std::vector<double>::iterator &it);
	std::vector<double>::iterator get_accumulators(std::vector<double>::iterator &it);
	std::vector<double>::iterator get_accumulatorRates(std::vector<double>::iterator &it);
	void print(std::ostream &out = std::cout);

};



#endif 
