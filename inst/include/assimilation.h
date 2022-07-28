#ifndef PLANT_FATE_PLANT_ASSIMILATION_H_
#define PLANT_FATE_PLANT_ASSIMILATION_H_

#include <phydro.h>

#include "plant_params.h"
#include "plant_geometry.h"

namespace plant{

struct PlantAssimilationResult{
	double gpp = 0;
	double npp = 0;
	double trans = 0;

	double dpsi_avg = 0;
	double vcmax_avg = 0;
	double gs_avg = 0;
	double c_open_avg = 0;

	double rleaf = 0;
	double rroot = 0;
	double rstem = 0;

	double tleaf = 0; 
	double troot = 0;
	
};

class Assimilator{
	public:
	// ~~ Last calculated values of rates
	// ~~ These are defined here rather than in local scope for debugging purposes. 
	PlantAssimilationResult plant_assim;
	// ~~

	public:	
	// **
	// ** Gross and Net Assimilation 
	// **
	template<class _Climate>
	phydro::PHydroResult leaf_assimilation_rate(double I0, double fapar, _Climate &clim, PlantParameters &par, PlantTraits &traits);
	
	template<class Env>
	void  calc_plant_assimilation_rate(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits);

	template<class Env>
	PlantAssimilationResult net_production(Env &env, PlantGeometry *G, PlantParameters &par, PlantTraits &traits);

	// ** 
	// ** Respiration and turnover
	// **
	//// leaf respiration rate - should be calculated AFTER asimialtion (needs updated Phydro outputs)
	double leaf_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits);
	double root_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits);
	double sapwood_respiration_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits);

	double leaf_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits);
	double root_turnover_rate(PlantGeometry *G, PlantParameters &par, PlantTraits &traits);
	

};

} // namespace plant

#include "assimilation.tpp"

#endif


