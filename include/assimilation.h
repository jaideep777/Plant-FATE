#ifndef PLANT_FATE_PLANT_ASSIMILATION_H_
#define PLANT_FATE_PLANT_ASSIMILATION_H_

#include <cmath>

#include "plant_params.h"
#include "plant_geometry.h"

namespace plant{

class Assimilator : public PlantGeometry{
	
	// **
	// ** Gross and Net Assimilation 
	// **
	template<class Env>
	double leaf_assimilation_rate(double fapar, Env &env){
		double P = 0.75;	// kg/m2/yr
		return P;
	}
	
	template<class Env>
	double plant_assimilation_rate(double fapar, Env &env){
		return leaf_assimilation_rate(fapar, env) * leaf_area;
	}

	template<class Env>
	double biomass_growth_rate(double fapar, Env &env, PlantParameters &par, PlantTraits &traits){
		double A = plant_assimilation_rate(fapar, env)*par.cbio;	// mol CO2 yr-1 * kg/mol CO2 = kg yr-1
		double R = 0.1*A // leaf_respiration_rate(G,par,traits)		// kg yr-1  [FIXME: For now, leaf respiration rate is simply 10% of production]
				 + root_respiration_rate(par,traits)		// kg yr-1
				 + sapwood_respiration_rate(par,traits);  // kg yr-1
		double T = leaf_turnover_rate(par,traits)			// kg yr-1
				 + root_turnover_rate(par,traits);		// kg yr-1

		return par.y*(A-R) - T;	// net biomass growth rate (kg yr-1)
	}

	// ** 
	// ** Respiration and turnover
	// ** 
	double leaf_respiration_rate(PlantParameters &par, PlantTraits &traits){
		return par.rl * traits.vcmax * leaf_area;
	}

	double root_respiration_rate(PlantParameters &par, PlantTraits &traits){
		return par.rr * root_mass(par, traits);
	}

	double sapwood_respiration_rate(PlantParameters &par, PlantTraits &traits){
		return par.rs * sapwood_mass(par, traits);
	}

	double leaf_turnover_rate(PlantParameters &par, PlantTraits &traits){
		return leaf_mass(par, traits) / traits.ll;	
	}
	
	double root_turnover_rate(PlantParameters &par, PlantTraits &traits){
		return root_mass(par, traits) / par.lr;
	}

	
};


} // namespace plant

#endif


