#include "assimilation.h"

#include <cmath>

namespace plant{


void Assimilator::les_update_lifespans(double lai, PlantParameters &par, PlantTraits &traits){
	double hT = plant_assim.vcmax_avg / plant_assim.vcmax25_avg;
	double f = 1;
	double fac = sqrt(((par.les_k1 * par.les_k2)*(par.les_k1 * par.les_k2) * f * hT * plant_assim.mc_avg) / (2 * par.les_u * par.les_cc));
	
	kappa_l = 365 * plant_assim.vcmax25_avg / (traits.lma*1e3) * fac;
	kappa_r = 365 * plant_assim.vcmax25_avg / (0.1333*1e3) * fac;
	//kappa_r = kappa_l * (par.les_cc/lai - 1) / (traits.zeta / traits.lma);
}


double Assimilator::les_assim_reduction_factor(phydro::PHydroResult& res, PlantParameters &par){
	double hT = res.vcmax / res.vcmax25;
	double f = 1;
	return 1; // Not applying age-related reduction factor because Phydro is already calibrated for average leaves (not yound leaves)
	// return 1 - sqrt(par.les_cc / (2 * par.les_u * res.mc * hT * f));
}


//// leaf respiration rate - should be calculated AFTER asimialtion (needs updated Phydro outputs)
double Assimilator::leaf_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//double vcmax_kg_yr = photo_leaf.vcmax * par.cbio * G->leaf_area;  // mol-CO2 m-2 year-1 * kg / mol-CO2 * m2
	//return par.rd * vcmax_kg_yr;
	return plant_assim.rleaf; // + par.rl * G->leaf_mass(traits);
}


double Assimilator::root_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return par.rr * G->root_mass(traits) * (plant_assim.gpp/G->crown_area);
}


double Assimilator::sapwood_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//return par.rs * G->sapwood_mass(traits);
//	double dpsi_gravity = (1000*10*G->height/1e6);
	double factor = traits.p50_xylem;
	double factor1 = (1 + 0.02*factor*factor); // 3e3 7e3
	return par.rs * G->sapwood_mass(traits)*factor1; // * (plant_assim.gpp/G->crown_area);
}


double Assimilator::leaf_turnover_rate(double _kappa_l, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return G->leaf_mass(traits) * _kappa_l; // / traits.ll;	
}


double Assimilator::root_turnover_rate(double _kappa_r, PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return G->root_mass(traits) * _kappa_r; // / par.lr;
}



} // namespace plant



