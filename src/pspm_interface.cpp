#include "pspm_interface.h"
#include "trait_evolution.h"
#include <iomanip>
#include "pspm_dynamic_environment.h"
using namespace std;

typedef PSPM_Dynamic_Environment EnvUsed;

// **********************************************************
// ********** PSPM_Plant ************************************
// **********************************************************

PSPM_Plant::PSPM_Plant() : plant::Plant() {
	
}

void PSPM_Plant::set_size(const std::array<double, STATE_DIM>& _x){
	Plant::set_size(_x[0]); // since Plant already has a set_size()! Keeping this overridden function for completeness. 
}

// This gives initial density of stages of the tree (established seedlings and thereafter)
double PSPM_Plant::init_density(void * _env, double input_seed_rain){
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
	return 1e-2*exp(-x[0]/0.1);
}


void PSPM_Plant::preCompute(double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	calc_demographic_rates(*env, t);
//	double p_plant_survival = exp(-vars.mortality);
//	//viable_seeds_dt = vars.fecundity_dt; // only for single-plant testrun
//	viable_seeds_dt = vars.fecundity_dt * p_plant_survival * env->patch_survival(t) / env->patch_survival(t_birth);
}

void PSPM_Plant::afterStep(double t, void * _env){
//	EnvUsed * env = (EnvUsed*)_env;
//	
//	seeds_hist.push(t, rates.dseeds_dt_germ);
	//seeds_hist.print_summary();
}

// Probability that a fresh seed survives to become a seedling
double PSPM_Plant::establishmentProbability(double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	return p_survival_dispersal(env) * p_survival_germination(*env);
}


std::array<double,STATE_DIM> PSPM_Plant::growthRate(double t, void * _env){
	return {rates.dsize_dt};
}

double PSPM_Plant::mortalityRate(double t, void * _env){
	double mort = rates.dmort_dt;
//	double mort_cndd = 
	return mort;
}

double PSPM_Plant::birthRate(double t, void * _env){
//	if (par.T_seed_rain_avg > 0){ 
//		return seeds_hist.get(); // birth rate is moving average of rate of germinating seeds over successional cycles
//	}          
//	else{
		return rates.dseeds_dt;       // fecundity rate is instantaneous rate of seed production 
//	}
}


// FIXME: This is used for initialization of every new cohort, not just the initial ones!
// So this must be consistent with init_density()
// Therefore needs size as input
void PSPM_Plant::init_accumulators(double t, void * _env){
	//set_size(x);	
	EnvUsed * env = (EnvUsed*)_env;
	geometry.lai = par.lai0;
	state.mortality = -log(establishmentProbability(t, env)); ///env->patch_survival(t));    // mortality // TODO: Verify! This is supposed to be cumulative mortality starting from the fresh seed stage until seedling stage
//	state.seed_pool = 0; // viable seeds
	t_birth = t;			// set cohort's birth time to current time
	//vars.mortality = 0; // only for single plant testrun
}

vector<double>::iterator PSPM_Plant::set_accumulators(vector<double>::iterator &it){
	geometry.lai       = *it++;
	state.mortality    = *it++;
//	state.seed_pool    = *it++;
	//vars.fecundity = viable_seeds; // only for single plant test run
	return it;
}

vector<double>::iterator PSPM_Plant::get_accumulators(vector<double>::iterator &it){
	*it++ = geometry.lai;
	*it++ = state.mortality;
//	*it++ = state.seed_pool;
	return it;
}

vector<double>::iterator PSPM_Plant::get_accumulatorRates(vector<double>::iterator &it){

	*it++ = rates.dlai_dt;	// lai
	*it++ = rates.dmort_dt; // mortality
//	*it++ = rates.dseeds_dt_pool; // seed pool size
	return it;
}


/// @brief      Print out evolvable traits and other important variables of the plant in a single line
/// @param out  stream to print to. 
/// @ingroup    trait_evolution 
/// @details    This function is called by the solver when printing a species. 
void PSPM_Plant::print(std::ostream &out) const {
	out << std::setw(10) << setprecision(3) << traits.species_name;
	vector<double> traits_vec = (const_cast<PSPM_Plant*>(this))->get_evolvableTraits();
	for (auto e : traits_vec){
		out << std::setw(10) << setprecision(5) << "|" << e << "| ";
	}
	out << std::setw(10) << setprecision(3) << geometry.get_size() 
	    << std::setw(10) << setprecision(3) << rates.dsize_dt 
	    << std::setw(10) << setprecision(3) << geometry.lai 
	    << std::setw(10) << setprecision(3) << state.mortality 
	    << std::setw(10) << setprecision(3) << rates.dseeds_dt 
	    << std::setw(10) << setprecision(3) << traits.a 
	    << std::setw(10) << setprecision(3) << traits.c 
	    ;
}


void PSPM_Plant::save(std::ostream &fout){
}

void PSPM_Plant::restore(std::istream &fin){
}

