#include "assimilation.h"

#include <cmath>

namespace plant{

// ** 
// ** Respiration and turnover
// **
//// leaf respiration rate - should be calculated AFTER asimialtion (needs updated Phydro outputs)
double Assimilator::leaf_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	//double vcmax_kg_yr = photo_leaf.vcmax * par.cbio * G->leaf_area;  // mol-CO2 m-2 year-1 * kg / mol-CO2 * m2
	//return par.rl * vcmax_kg_yr;
	return plant_assim.rleaf;
}

double Assimilator::root_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return par.rr * G->root_mass(traits);
}

double Assimilator::sapwood_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return par.rs * G->sapwood_mass(traits);
}

double Assimilator::leaf_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return G->leaf_mass(traits) / traits.ll;	
}

double Assimilator::root_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits){
	return G->root_mass(traits) / par.lr;
}



} // namespace plant



