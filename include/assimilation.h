#ifndef PLANT_FATE_PLANT_ASSIMILATION_H_
#define PLANT_FATE_PLANT_ASSIMILATION_H_

#include <cmath>

#include "plant_params.h"
#include "plant_geometry.h"

namespace plant{

class Assimilator{
	public:
	// ~~ Last calculated values of rates
	// ~~ These are defined here rather than in local scope for debugging purposes. 
	double A,R,T, rl, rr, rs, tl, tr, Anet;
	// ~~
	
	public:	
	// **
	// ** Gross and Net Assimilation 
	// **
	template<class Env>
	double leaf_assimilation_rate(double fapar, Env &env){
		double P = 15;	// umol/m2/s
		return P * (1e-6*86400*365.2524);	// umol m-2 s-1 --> mol m-2 yr-1
	}
	
	template<class Env>
	double plant_assimilation_rate(double fapar, Env &env, PlantGeometry *G){
		return leaf_assimilation_rate(fapar, env) * G->leaf_area;
	}

	template<class Env>
	double biomass_growth_rate(double fapar, Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		A = plant_assimilation_rate(fapar, env, G) * par.cbio;	// mol CO2 yr-1 * kg/mol CO2 = kg yr-1
		
		rl = leaf_respiration_rate(G,par,traits);		// kg yr-1  [FIXME: For now, leaf respiration rate is simply 10% of production]
		rr = root_respiration_rate(G, par,traits);		// kg yr-1
		rs = sapwood_respiration_rate(G, par,traits);  // kg yr-1
		
		tl = leaf_turnover_rate(G, par,traits);			// kg yr-1
		tr = root_turnover_rate(G, par,traits);		// kg yr-1
		
		R = rl+rr+rs;
		T = tl+tr;

		Anet = par.y*(A-R) - T;	// net biomass growth rate (kg yr-1)
		return Anet;
	}

	// ** 
	// ** Respiration and turnover
	// **
	// leaf respiration rate - should be calculated AFTER asimialtion (needs updated Vcmax)
	double leaf_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		double vcmax_kg_yr = traits.vcmax * par.cbio * G->leaf_area;  // mol-CO2 m-2 year-1 * kg / mol-CO2 * m2
		return par.rl * vcmax_kg_yr;
	}

	double root_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		return par.rr * G->root_mass(par, traits);
	}

	double sapwood_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		return par.rs * G->sapwood_mass(par, traits);
	}

	double leaf_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		return G->leaf_mass(par, traits) / traits.ll;	
	}
	
	double root_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		return G->root_mass(par, traits) / par.lr;
	}


};


} // namespace plant

#endif


