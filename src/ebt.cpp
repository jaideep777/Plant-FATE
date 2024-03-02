#include <algorithm>
#include <cassert>

#include <solver.h>
using namespace std;


void Solver::calcRates_EBT(double t, vector<double>::iterator S, vector<double>::iterator dSdt){

	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;

	for (int s = 0; s < species_vec.size(); ++s){
		Species_Base * spp = species_vec[s];

		std::vector<double> pi0 = spp->getX(spp->J-1);	 // last cohort is pi0, N0
		double N0 = spp->getU(spp->J-1);
		//std::cout << "pi = " << pi0 << ", N0 = " << N0 << "\n";

		std::vector<double> grad_dx(spp->istate_size, control.ebt_grad_dx);

		std::vector<std::vector<double>> g_gx = spp->growthRateGradient(-1, t, env, grad_dx);
		std::vector<double> m_mx = spp->mortalityRateGradient(-1, t, env, grad_dx);
		//std::cout << "g = " << g_gx[0] << ", gx = " << g_gx[1] << "\n";
		//std::cout << "m = " << m_mx[0] << ", mx = " << m_mx[1] << "\n";

		double birthFlux;
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * pe;
		}
		else{
			birthFlux = spp->birth_flux_in * pe;
		}

		for (int i=0; i<spp->J-1; ++i){	// go down to the second last cohort (exclude boundary cohort)
			std::vector<double> dx = spp->growthRate(i, t, env);							// dx/dt
			double du = -spp->mortalityRate(i, t, env) * spp->getU(i);		// du/dt
			//std::cout << "S/C = " << s << "/" << i << " " << spp->getX(i) << " " << dx << " " << du << "\n";
			// add all the rates to the iterator
			for(int k = 0; k < dx.size(); ++k) *itr++ = dx[k];
			*itr++ = du;
			// not sure what this is actually - JJ: Just keeping the rate and state iterators at equivalent positions, though it might not be necassary
			its += 1 + dx.size(); 
		}

		// dpi0/dt and dN0/dt
		//std::cout << "S/C = " << s << "/" << "b" << " | pi0/N0 = " << pi0 << " " << N0;
		//if (pi0 <= 0) pi0 = 1e-40;
		//if (N0  <= 0) N0  = 1e-40;
		vector<double> dpi0(spp->istate_size);
		for (int k=0; k<spp->istate_size; ++k){
			dpi0[k] = g_gx[0][k]*N0 - m_mx[0]*pi0[k];
			for (int j=0; j<spp->istate_size; ++j){
				dpi0[k] += g_gx[j+1][k]*pi0[j];
			}
			*itr++ = dpi0[k];
			its++;
			//std::cout << " | dpi0/dN0 = " << dpi0 << " " << dN0 << " | mx/gx/mb/gb = " << m_mx[1] << " " << g_gx[1] << " " << m_mx[0] << " " << g_gx[0] << "\n";
		}

		double dN0  = -m_mx[0]*N0 - std::inner_product(pi0.begin(), pi0.end(), std::next(m_mx.begin()), 0.0) + birthFlux;
		*itr++ = dN0;
		its++;

		if (spp->n_accumulators > 0){
			auto itr_prev = itr;
			spp->accumulatorRates(itr);
			assert(distance(itr_prev, itr) == spp->n_accumulators*spp->J); 
			its += spp->n_accumulators*spp->J; 	
		}
		
	}

}



void Solver::addCohort_EBT(){
	// Q: is updateEnv needed here (and in cm.cpp) to init cumulative vars?
	// Maybe not, see logic in abm.cpp
	for (auto spp : species_vec){
		// 1. internalize the pi0-cohort (this cohort's birth time would have been already set when it was inserted)
		realizeEbtBoundaryCohort(spp);

		// 2. insert a new cohort (copy of boundary cohort, to be the new pi0-cohort)
		spp->initBoundaryCohort(current_time, env);	// update initial extra state and birth-time of boundary cohort
		spp->addCohort(); // introduce copy of boundary cohort into species

		// 3. set x,u of the new cohort to 0,0, thus marking it as the new pi0-cohort
		spp->setX(spp->J-1, std::vector<double>(spp->istate_size, 0)); 
		spp->setU(spp->J-1, 0);
	}

	// 4. reset state from cohorts
	resizeStateFromSpecies(); 
	copyCohortsToState();
	
}


void Solver::removeDeadCohorts_EBT(){

	for (auto spp : species_vec){
		spp->markDeadCohorts(control.ebt_ucut);
		spp->removeMarkedCohorts();
	}

	resizeStateFromSpecies();
	copyCohortsToState();

}

void Solver::mergeCohorts_EBT(){
	// for (auto spp : species_vec){
	// 	spp->mergeCohortsAddU(control.ebt_merge_dxcut);
	// }

	// resizeStateFromSpecies();
	// copyCohortsToState();

}
