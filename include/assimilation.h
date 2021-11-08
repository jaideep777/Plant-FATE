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
	double A,R,T, rl, rr, rs, tl, tr, Anet, rl_temp;
	// ~~
	//phydro::PHydroResult photo_leaf;
	//double E_per_leaf_area;
	//double dpsi_x;
	//double psi_x;

	public:	
	// **
	// ** Gross and Net Assimilation 
	// **
	template<class _Climate>
	phydro::PHydroResult leaf_assimilation_rate(double fapar, _Climate &clim, PlantParameters &par, PlantTraits &traits){
		phydro::ParCost par_cost(par.alpha, par.gamma);
		phydro::ParPlant par_plant(traits.K_leaf, traits.p50_leaf, traits.b_leaf);
		par_plant.gs_method = phydro::GS_APX;
		auto photo_leaf = phydro::phydro_analytical(clim.tc,  clim.ppfd*4,  clim.vpd,  clim.co2,	// FIXME: ppfd*4 here to convert daily average PAR to daily max PAR.
													clim.elv,    fapar,    par.kphio,  clim.swp, 
													par.rl, par_plant, par_cost);
		
		return photo_leaf;	// umol m-2 s-1 
	}
	
	template<class Env>
	double plant_assimilation_rate(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		double GPP_plant = 0, Rl_plant = 0;
		
		double lai = par.lai_max*traits.fl;
		double fapar = 1-exp(-0.5*lai);
		
		auto res = leaf_assimilation_rate(fapar, env.clim, par, traits);
		double la_layer = G->leaf_area;
		
		GPP_plant += (res.a + res.vcmax*par.rl) * la_layer;
		Rl_plant  += (res.vcmax*par.rl) * la_layer;
	
		// calculate yearly averages in mol/yr	
		double f_light_day = 0.25; // fraction day that receives max light (x0.5 sunlight hours, x0.5 average over sinusoid)
		double f_growth_yr = 1.0;  // factor to convert daily mean PAR to yearly mean PAR
		double f = f_light_day * f_growth_yr * (86400*365.2524) * 1e-6; // umol/s ---> mol/yr

		rl_temp = Rl_plant * f;	
		return GPP_plant * f;
	}

	template<class Env>
	double biomass_growth_rate(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		A = plant_assimilation_rate(env, G, par, traits) * par.cbio;	// mol CO2 yr-1 * kg/mol CO2 = kg yr-1
		
		rl = leaf_respiration_rate(G,par,traits);		// kg yr-1  
		rr = root_respiration_rate(G, par,traits);		// kg yr-1
		rs = sapwood_respiration_rate(G, par,traits);	// kg yr-1
		
		tl = leaf_turnover_rate(G, par,traits);			// kg yr-1
		tr = root_turnover_rate(G, par,traits);		// kg yr-1
		
		R = rl + rr + rs;
		T = tl + tr;

		Anet = par.y*(A-R) - T;	// net biomass growth rate (kg yr-1)
		return Anet;
	}

	// ** 
	// ** Respiration and turnover
	// **
	//// leaf respiration rate - should be calculated AFTER asimialtion (needs updated Phydro outputs)
	double leaf_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
		//double vcmax_kg_yr = photo_leaf.vcmax * par.cbio * G->leaf_area;  // mol-CO2 m-2 year-1 * kg / mol-CO2 * m2
		//return par.rl * vcmax_kg_yr;
		return rl_temp * par.cbio;
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


