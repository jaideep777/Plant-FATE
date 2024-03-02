#include "pspm_plant.h"
using namespace std;

typedef LightEnvironment EnvUsed;


PSPM_Plant::PSPM_Plant() : plant::Plant() {
	
}

void PSPM_Plant::set_size(const std::array <double, 1>& _x){
	set_height(_x[0]);
}

double PSPM_Plant::init_density(void * _env, double input_seed_rain){
	//if (x == seed.vars.height){
		//p.set_height(x);
		EnvUsed * env = (EnvUsed*)_env;
		compute_vars_phys(*env);
		double u0;
#ifndef USE_INIT_DIST
		if (x[0] < 0.7)
			u0 = fabs(input_seed_rain)*germination_probability(*env)/vars.height_dt;
		else 
			u0 = 0;
		return u0;
#else
		// if (lma == 0.0825) u0 = 1e-2;
		// else if (lma == 0.02625) u0 = 1e-3;
		// else if (lma == 0.04625) u0 = 1e-4;
		u0 = 1e-2*germination_probability(*env)/vars.height_dt;
		return (x[0]>15)? 1e-16 : u0*exp(-5*x[0]/20);
#endif
	//}
	//else return 0;
}


void PSPM_Plant::preCompute(double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	compute_vars_phys(*env);
	double p_plant_survival = exp(-vars.mortality);
	//viable_seeds_dt = vars.fecundity_dt; // only for single-plant testrun
	viable_seeds_dt = vars.fecundity_dt * p_plant_survival * env->patch_survival(t) / env->patch_survival(t_birth);
}

double PSPM_Plant::establishmentProbability(double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	return germination_probability(*env);
}

std::array<double,1> PSPM_Plant::growthRate(double t, void * env){
	//if (p.vars.height != x){
		//p.set_height(x);
		//compute_vars_phys(*env);
		//++nrc;
	////}
	//cout << "x/g = " << x << " " << vars.height_dt << endl;
	return {vars.height_dt};
		
}

double PSPM_Plant::mortalityRate(double t, void * env){
	//assert(p.vars.height == x);
	//if (p.vars.height != x){
		//p.set_height(x);
		//p.compute_vars_phys(*env);
		//++ndc;
	//}
	return vars.mortality_dt;
}

double PSPM_Plant::birthRate(double t, void * env){
	// Need this here because birthRate is not called in order, and only called rarely,
	// after completion of step_to.
	//if (p.vars.height != x){
		//++nbc;
		//p.set_height(x);
		//p.compute_vars_phys(*env);
	//}
	//assert(p.vars.height == x);
	return vars.fecundity_dt;
}


void PSPM_Plant::init_accumulators(double t, void * _env){
	//set_size(x);	
	EnvUsed * env = (EnvUsed*)_env;
	vars.mortality = -log(germination_probability(*env)); ///env->patch_survival(t));    // mortality
	viable_seeds = 0; // viable seeds
	vars.area_heartwood = 0;   // heartwood area
	vars.mass_heartwood = 0;     // heartwood mass
	t_birth = t;			// set cohort's birth time to current time
	//vars.mortality = 0; // only for single plant testrun
}

vector<double>::iterator PSPM_Plant::set_accumulators(vector<double>::iterator &it){
	vars.mortality      = *it++;
	viable_seeds        = *it++;
	vars.area_heartwood = *it++;
	vars.mass_heartwood = *it++;
	//vars.fecundity = viable_seeds; // only for single plant test run
	return it;
}

vector<double>::iterator PSPM_Plant::get_accumulators(vector<double>::iterator &it){
	*it++ = vars.mortality;
	*it++ = viable_seeds;
	*it++ = vars.area_heartwood;
	*it++ = vars.mass_heartwood;
	return it;
}

vector<double>::iterator PSPM_Plant::get_accumulatorRates(vector<double>::iterator &it){

	*it++ = vars.mortality_dt;	// mortality
	*it++ = viable_seeds_dt; // viable_seeds
	*it++ = vars.area_heartwood_dt; // heartwood area
	*it++ = vars.mass_heartwood_dt; // heartwood mass
	return it;
}

void PSPM_Plant::print(std::ostream &out){
	out << "|" << lma << "|\t" << vars.mortality << "\t" << vars.fecundity << "\t" << viable_seeds << "\t" << vars.area_heartwood << "\t" << vars.mass_heartwood << "\t";
}

