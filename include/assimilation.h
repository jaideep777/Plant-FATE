#ifndef PLANT_FATE_PLANT_ASSIMILATION_H_
#define PLANT_FATE_PLANT_ASSIMILATION_H_

#include <cmath>

#include <phydro.h>

#include "plant_params.h"
#include "plant_geometry.h"

namespace plant{

class Assimilator{
	public:
	// ~~ Last calculated values of rates
	// ~~ These are defined here rather than in local scope for debugging purposes. 
	double A,R,T, rl, rr, rs, tl, tr, Anet;
	// ~~
	phydro::PHydroResult photo_leaf;
	double E_per_leaf_area;
	double dpsi_x;
	double psi_x;

	public:	
	// **
	// ** Gross and Net Assimilation 
	// **
	template<class Env>
	double leaf_assimilation_rate(double fapar, Env &env, PlantParameters &par, PlantTraits &traits){
		phydro::ParCost par_cost(par.alpha, par.gamma);
		phydro::ParPlant par_plant(traits.K_leaf, traits.p50_leaf, traits.b_leaf);
		par_plant.gs_method = phydro::GS_APX;
		photo_leaf = phydro::phydro_analytical(env.tc,  env.ppfd,  env.vpd,   env.co2, 
											   env.elv,    fapar,  par.kphio, env.swp, 
											   par.rl, par_plant, par_cost);
		double A_inst = photo_leaf.a;
		//std::cout << "A_inst = " << A_inst << "\n";
		//double A_inst = 40;
	   	double A_avg_day = A_inst * 86400;
		double A_avg_yr = A_avg_day * 365.2524;	// half the year is growth season 
		return A_avg_yr * (1e-6);	// umol m-2 yr-1 --> mol m-2 yr-1
	}
	
	template<class Env>
	double plant_assimilation_rate(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		return leaf_assimilation_rate(1.0, env, par, traits) * G->leaf_area;
	}

	template<class Env>
	double biomass_growth_rate(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		A = plant_assimilation_rate(env, G, par, traits) * par.cbio;	// mol CO2 yr-1 * kg/mol CO2 = kg yr-1
		
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


