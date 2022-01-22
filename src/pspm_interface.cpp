#include "pspm_interface.h"
using namespace std;

typedef PSPM_Dynamic_Environment EnvUsed;


PSPM_Plant::PSPM_Plant() : plant::Plant() {
	
}

void PSPM_Plant::set_size(double _x){
	Plant::set_size(_x); // since Plant already has a set_size()! Keeping this overridden function for completeness. 
}

double PSPM_Plant::init_density(double x, void * _env, double input_seed_rain){
//	EnvUsed * env = (EnvUsed*)_env;
//	compute_vars_phys(*env);
//	double u0;
//	if (x == vars.height){
//		if (input_seed_rain < 0) input_seed_rain = 1;
//		cout << "Here, setting u0 (" << input_seed_rain << "). " << vars.height << " " << x << "\n";
//		u0 = input_seed_rain*germination_probability(*env)/vars.height_dt;
//	}
//	else{
//		cout << "Here, setting u0 to 0. " << vars.height << " " << x << "\n";
//		u0 = 0;
//	}
//	return u0;
//	if (geometry.diameter < 0.02) return 1;
//	else return 0; // dummy initial density of 1, for now. This shouldn't matter because there is spinup.
	return 0;
}


void PSPM_Plant::preCompute(double x, double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	calc_demographic_rates(*env);
//	double p_plant_survival = exp(-vars.mortality);
//	//viable_seeds_dt = vars.fecundity_dt; // only for single-plant testrun
//	viable_seeds_dt = vars.fecundity_dt * p_plant_survival * env->patch_survival(t) / env->patch_survival(t_birth);
}

void PSPM_Plant::afterStep(double x, double t, void * _env){
//	EnvUsed * env = (EnvUsed*)_env;
//	
//	seeds_hist.push(t, rates.dseeds_dt_germ);
	//seeds_hist.print_summary();
}


double PSPM_Plant::establishmentProbability(double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	return p_survival_germination(*env);
}


double PSPM_Plant::growthRate(double x, double t, void * _env){
	return rates.dsize_dt;
}

double PSPM_Plant::mortalityRate(double x, double t, void * _env){
	return rates.dmort_dt;
}

double PSPM_Plant::birthRate(double x, double t, void * _env){
//	if (par.T_seed_rain_avg > 0){ 
//		return seeds_hist.get(); // birth rate is moving average of rate of germinating seeds over successional cycles
//	}          
//	else{
		return rates.dseeds_dt_germ;       // birth rate is instantaneous rate of germinating seeds 
//	}
}

void PSPM_Plant::init_state(double t, void * _env){
	//set_size(x);	
	EnvUsed * env = (EnvUsed*)_env;
	geometry.lai = par.lai0;
	geometry.crootmass = 0;
	state.mortality = -log(p_survival_germination(*env)); ///env->patch_survival(t));    // mortality
	state.seed_pool = 0; // viable seeds
	t_birth = t;			// set cohort's birth time to current time
	//vars.mortality = 0; // only for single plant testrun
}

vector<double>::iterator PSPM_Plant::set_state(vector<double>::iterator &it){
	geometry.lai    = *it++;
	state.mortality  = *it++;
	state.seed_pool  = *it++;
	//vars.fecundity = viable_seeds; // only for single plant test run
	return it;
}

vector<double>::iterator PSPM_Plant::get_state(vector<double>::iterator &it){
	*it++ = geometry.lai;
	*it++ = state.mortality;
	*it++ = state.seed_pool;
	return it;
}

vector<double>::iterator PSPM_Plant::get_rates(vector<double>::iterator &it){

	*it++ = rates.dlai_dt;	// lai
	*it++ = rates.dmort_dt; // mortality
	*it++ = rates.dseeds_dt_pool; // seed pool size
	return it;
}

void PSPM_Plant::print(std::ostream &out){
	out << traits.species_name << "\t"
		<< "|" << traits.lma << "|\t" 
	    << geometry.get_size() << "\t" 
	    << rates.dsize_dt << "\t" 
	    << geometry.lai << "\t" 
	    << state.mortality << "\t" 
	    << state.seed_pool << "\t" 
	    ;
}





