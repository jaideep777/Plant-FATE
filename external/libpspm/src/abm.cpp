#include <cassert>
#include "solver.h"
using namespace std;

void Solver::stepABM(double t, double dt){
	if (debug) std::cout << "Entering stepABM ..." << std::endl;
	// Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		if (debug) std::cout << "spp " << s << ": J = " << spp->J << std::endl;
		// Calc growth and mortality rates
		vector<vector<double>> growthArray(spp->J);
		vector<double> mortalityArray(spp->J);
		for (int i=0; i<spp->J; ++i){
			growthArray[i]    = spp->growthRate(i, t, env);
			mortalityArray[i] = spp->mortalityRate(i, t, env);
		}

		// calculate fecundity
		double birthFlux;
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * pe;
		}
		else{
			birthFlux = spp->birth_flux_in * pe;
		}
		
		double noff = birthFlux*dt/spp->get_boundary_u();  // number of offspring = birthflux / density of each superindividual
		if (debug) std::cout << "noff: " << spp->noff_abm << " | " << noff << " " << birthFlux << " " << dt << " " << spp->get_boundary_u() << std::endl;
		spp->noff_abm += noff;
		
		// implement growth
		for (int i=0; i<spp->J; ++i){
			std::vector<double> X = spp->getX(i);
			for(int k = 0; k<X.size(); ++k){
				X[k] += growthArray[i][k] * dt;
			}
			spp->setX(i, X);
		}

		// implement mortality
		for (int i=0; i<spp->J; ++i){
			double survival_prob = exp(-mortalityArray[i]*dt);
			//cout << "survival prob = " << survival_prob << "\n";
			if (double(rand())/RAND_MAX >= survival_prob) spp->markCohortForRemoval(i);
		}
		spp->removeMarkedCohorts();
	}

	// step extra variables
	// these will be in order abc abc abc... 
	// use the state vector itself to temporarily store and retrive state, since it's not used for x and u (but after system variables)

	// Elisa: copied over from ebt.cpp after code reformatting
	
	for (auto spp : species_vec){
		vector<double>::iterator its = state.begin() + n_statevars_system; // Skip system variables
		vector<double>::iterator itr = rates.begin() + n_statevars_system;
		if (spp->n_accumulators > 0){
			auto itr_prev = itr;
			auto its_prev = its;

			spp->copyAccumulatorsToState(its); // this will increment its
			spp->accumulatorRates(itr);      // this will increment itr
			assert(distance(itr_prev, itr) == spp->n_accumulators*spp->J);
			assert(distance(its_prev, its) == spp->n_accumulators*spp->J);

			itr = itr_prev; its = its_prev; // so bring them back
			
			for (int i=0; i< spp->n_accumulators*spp->J; ++i){
				*its += (*itr)*dt;
				++its; ++itr;
			}
			itr = itr_prev; its = its_prev; // so bring them back

			spp->copyAccumulatorsToCohorts(its); // this will increment its
		}
	}

	// Q: Should updated environment be used for initializing cummulative variables in new cohorts?
	// No, Should not be done, because recruitment is happening continuously, not necessarily after mortality.
	// It's best to do recruit initialization at the beginning of the timestep (i.e., with original env and at time t, istead of updating env here and using t+dt)
	// Also, using updateEnv here gives wrong results! 
	// updateEnv(t+dt, state.begin(), rates.begin()); // <-- don't do this

	// Insert offspring
	for (auto spp : species_vec){
		// add recruits once they have accumulated to > 1
		// use t and env to initialize (instead of updated values) becaues number of recruits were calculated at the beginning of the step
		if (debug) std::cout << "spp " << spp << ": adding cohorts; n = " << int(spp->noff_abm) << std::endl;

		if (spp->noff_abm > 1){
			int nadd = int(spp->noff_abm);	
			spp->initBoundaryCohort(t, env);
			spp->addCohort(nadd);
			spp->noff_abm -= nadd;
		}
		// also add a recruit in anticipation if the whole populaiton is dead
		if (spp->J == 0){
			spp->initBoundaryCohort(t, env);
			spp->addCohort(1);
			//--spp->noff_abm;
		}
		
	}
	
	// resize state - this matters only if extra istate is specified
	resizeStateFromSpecies();
	copyCohortsToState();
	//print();
	if (debug) std::cout << "Exiting stepABM ..." << std::endl;

}


