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
	return par.rr * G->root_mass(traits);
}

double Assimilator::sapwood_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//return par.rs * G->sapwood_mass(traits);
//	double dpsi_gravity = (1000*10*G->height/1e6);
	double factor = traits.p50_xylem;
	double factor1 = (1 + 0.02*factor*factor); // 3e3 7e3
	return par.rs * G->sapwood_mass(traits)*factor1;	
}

double Assimilator::leaf_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	// double hT = 1; //exp(par.les_c - par.les_hT_dH/par.les_molar_R/)
	// double fac = par.les_k1 * par.les_k2 * 1 * plant_assim.mc / (2 * par.les_u * les_cc);
	// double kappa = plant_assim.vcmax_avg / traits.lma;
	return G->leaf_mass(traits) / traits.ll;	
}

double Assimilator::root_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return G->root_mass(traits) / par.lr;
}



} // namespace plant



