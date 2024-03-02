#include <cassert>
#include <string>
#include "solver.h"
using namespace std;

void Solver::stepU_iCM(double t, vector<double> &S, vector<double> &dSdt, double dt){
		
	vector<double>::iterator its = S.begin() + n_statevars_system; // Skip system variables
	
	// 1. Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		// [S S S x u x u x u x u x u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *XU = &(*its); 
		int J = spp->J;

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
		// cout << "B = " << birthFlux << ", ub = " << birthFlux/(gb[0]+1e-20) << '\n';

		for (int i=0; i<spp->J; ++i){
			vector<vector<double>> g_gx = spp->growthRateGradient(i, t, env, control.cm_grad_dx);  

			for (int k=0; k<spp->istate_size; ++k){
				*its++ += g_gx[0][k]*dt;
			}

			double du = spp->mortalityRate(i, t, env); 
			for (int k=0; k<spp->istate_size; ++k){
				du += g_gx[k+1][k]; 
			}
			if (!control.cm_use_log_densities) *its++ /= 1 + du*dt;
			else *its++ -= du*dt;
		}

		its += J*spp->n_accumulators;

		// // all cohorts
		// for (int i=0; i<J; ++i){
		// 	XU[2*i+0] += spp->growthRate(i, t, env)[0]*dt;

		// 	std::vector<std::vector<double>> g_gx = spp->growthRateGradient(i, t, env, control.cm_grad_dx);
		// 	XU[2*i+1] /= 1 + g_gx[1][0]*dt + spp->mortalityRate(i, t, env)*dt;
		// }
	
		// its += J*(2+spp->n_accumulators);

	}
	
	assert(distance(S.begin(), its)==S.size());
}


