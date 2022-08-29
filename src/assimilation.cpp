#include "assimilation.h"

#include <cmath>

namespace plant{

// ** 
// ** Respiration and turnover
// **
//// leaf respiration rate - should be calculated AFTER asimialtion (needs updated Phydro outputs)
double Assimilator::leaf_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//double vcmax_kg_yr = photo_leaf.vcmax * par.cbio * G->leaf_area;  // mol-CO2 m-2 year-1 * kg / mol-CO2 * m2
	//return par.rd * vcmax_kg_yr;
	return plant_assim.rleaf; // + par.rl * G->leaf_mass(traits);
}

double Assimilator::root_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return par.rr * G->root_mass(traits) * (plant_assim.gpp/G->crown_area/4.5);
}

double Assimilator::sapwood_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//return par.rs * G->sapwood_mass(traits);
//	double dpsi_gravity = (1000*10*G->height/1e6);
	double factor = traits.p50_xylem;
	double factor1 = (1 + 0.02*factor*factor); // 3e3 7e3
	return par.rs * G->sapwood_mass(traits)*factor1 * (plant_assim.gpp/G->crown_area/4.5);	
}

double Assimilator::leaf_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	double hT = plant_assim.vcmax_avg / plant_assim.vcmax25_avg;
	double f = 1;
	double fac = sqrt(par.les_k1 * par.les_k2 * f * hT * plant_assim.mc_avg / (2 * par.les_u * par.les_cc));
	double kappa = plant_assim.vcmax25_avg / (traits.lma*1e3) * fac;
	// std::cout << "vcmax / vcmax25 / ll = " << plant_assim.vcmax_avg << " / " << plant_assim.vcmax25_avg << " / " << 1.0/(kappa*365) << std::endl;
	return G->leaf_mass(traits) * kappa * 365; // / traits.ll;	
}

double Assimilator::root_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	double hT = plant_assim.vcmax_avg / plant_assim.vcmax25_avg;
	double f = 1;
	double fac = sqrt(par.les_k1 * par.les_k2 * f * hT * plant_assim.mc_avg / (2 * par.les_u * par.les_cc));
	double kappa = plant_assim.vcmax25_avg / (traits.lma*1e3) * fac;
	return G->root_mass(traits) * kappa * 365; // / par.lr;
}



} // namespace plant



