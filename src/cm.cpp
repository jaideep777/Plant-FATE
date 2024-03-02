#include <algorithm>
#include <cassert>

#include "solver.h"
//#include "vector_insert.h"
using namespace std;

// state must be copied to cohorts before calling this function
void Solver::calcRates_CM(double t, vector<double>::iterator S, vector<double>::iterator dSdt){

	//auto ss = species_vec[0];	
	//cout << "svec: " << t << " " << S[0] << " " << S[1] << " " << S[2] << " " << S[3] << "\n";
	//cout << "coho: " << t << " " << ss->getX(0) << " " << ss->getU(0) << " " << ss->getX(1) << " " << ss->getU(1) << "\n";
	
	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;

	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		//cout << "calcRates(species = " << s << ")\n";
		
		// Boundary u is not used in rate calcs per se, but needed in size integral. Hence update
		// FIXME: This boundary condition works only for 1D state
		double birthFlux;
		vector<double> gb = spp->growthRate(-1, t, env);
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * pe;
		}
		else{
			birthFlux = spp->birth_flux_in * pe;
		}
		spp->set_ub(birthFlux/(gb[0]+1e-20));

		for (int i=0; i<spp->J; ++i){
			vector<vector<double>> g_gx = spp->growthRateGradient(i, t, env, control.cm_grad_dx);  
			for (int k=0; k<spp->istate_size; ++k){
				double dx =  g_gx[0][k];
				*itr++ = dx;
			}

			double du = -spp->mortalityRate(i, t, env); 
			for (int k=0; k<spp->istate_size; ++k){
				du += -g_gx[k+1][k]; 
			}
			if (!control.cm_use_log_densities) du *= spp->getU(i);
			*itr++ = du;	
			
			its += (spp->istate_size+1);
		}

		
		// extra rates
		if (spp->n_accumulators > 0){
			auto itr_prev = itr;
			spp->accumulatorRates(itr);
			assert(distance(itr_prev, itr) == spp->n_accumulators*spp->J); 
			its += spp->n_accumulators*spp->J; 	
		}
		//cout << "---\n";
	}
}


void Solver::addCohort_CM(){
	//cout << ".......... add cohorts ...............\n";
	updateEnv(current_time, state.begin(), rates.begin());
	
	for (int s = 0; s< species_vec.size(); ++s){
		auto spp = species_vec[s];
		
		// calculate u for boundary cohort
		vector<double> gb = spp->growthRate(-1, current_time, env);
		double pe = spp->establishmentProbability(current_time, env);
		double birthFlux;
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,current_time) * pe;
		}
		else{
			birthFlux = spp->birth_flux_in * pe;
		}
		// cout << "Adding cohort: gb = " << gb[0] << '\n';
		spp->initBoundaryCohort(current_time, env); // init extra state variables and birth time of the boundary cohort
		spp->set_ub(birthFlux/(gb[0]+1e-20)); // FIXME: This line works only for 1D state
		spp->addCohort();	// introduce copy of boundary cohort in system
	}
	
	resizeStateFromSpecies();
	copyCohortsToState();
}


void Solver::removeCohort_CM(){

	for (auto& spp : species_vec){
		if (control.cm_remove_cohorts){ // remove a cohort if number of cohorts in the species exceeds a threshold
			spp->markDeadCohorts(control.ebt_ucut);
			spp->markDenseCohorts(control.cm_dxcut);
			spp->removeMarkedCohorts();
		}
	}
	
	resizeStateFromSpecies();
	copyCohortsToState();
	//cout << "........................................\n";
	
}



//template<class Model, class Environment>
//double Solver<Model,Environment>::calc_u0_CM(){
	////// function to iterate
	////auto f = [this](double utry){
	////    // set u0 to given (trial) value
	////    state[xsize()] = utry;
	////    // recompute environment based on new state
	////    mod->computeEnv(current_time, state, this);
	////    // calculate birthflux by trapezoidal integration (under new environment)
	////    double birthFlux = integrate_x([this](double z, double t){return mod->birthRate(z,t);}, current_time, state, 1);
		
	////    double unext = birthFlux/mod->growthRate(xb, current_time);
	////    return unext;
	////};

	////double u0 = state[xsize()+1]; // initialize with u0 = u1
	////// iterate
	////double err = 100;
	////while(err > 1e-6){
	////    double u1 = f(u0);
	////    err = abs(u1 - u0);
	////    u0 = u1;
	////}
	////state[xsize()] = u0;
	//////cout << "u0 = " << u0 << endl;	
	////return state[xsize()];
	//return 0;
//}

