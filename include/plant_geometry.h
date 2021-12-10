#ifndef PLANT_FATE_PLANT_GEOMETRY_H_
#define PLANT_FATE_PLANT_GEOMETRY_H_

#include <cmath>

#include "plant_params.h"

namespace plant{

class PlantGeometry{
	private:
	struct{
		// geometry traits
		double m, n;    // crown shape paramaters
		double a;       // height-diameter allometry
		double c;       // crown area allometry
		double fg;      // upper canopy gap fraction
		
		// Precomputed Geometric parameters
		double eta_c;
		double pic_4a;
		double zm_H;
		double qm;

		// allocation
		double dmat;    // diameter at reproductive maturity
	} geom;

	public:
	// current state
	double lai;                      // leaf area index 
	double diameter;                     // basal diameter

	// variables calculated from state variables
	double height;                       // height
	double crown_area;                   // crown area
	double sapwood_fraction;             // fraction of stem that is sapwood
	double functional_xylem_fraction;    // fraction of funcitonal xylem in sapwood

	// ode-based calculations of sapwood and heartwood (for debug)
	double sap_frac_ode = 1;
	double sapwood_mass_ode = 0;
	double heart_mass_ode = 0;
	double k_sap;

	public:

	void init(PlantParameters &par, PlantTraits &traits);

	// **
	// ** Crown geometry
	// **
	double q(double z);
	double zm();

	double crown_area_extent_projected(double z, PlantTraits &traits);
	double crown_area_above(double z, PlantTraits &traits);

	
	// **
	// ** Biomass partitioning
	// **
	double dsize_dmass(PlantTraits &traits) const ;

	double dreproduction_dmass(PlantParameters &par, PlantTraits &traits);


	// **
	// ** LAI model
	// ** 
	double dmass_dt_lai(double &dL_dt, double dmass_dt_max, PlantTraits &traits);


	// **
	// ** Carbon pools
	// **
	double leaf_mass(PlantTraits &traits);
	double root_mass(PlantTraits &traits);
	double sapwood_mass(PlantTraits &traits);
	double sapwood_mass_real(PlantTraits &traits);
	double stem_mass(PlantTraits &traits);
	double heartwood_mass(PlantTraits &traits);
	double total_mass(PlantTraits &traits);

	// **
	// ** state manipulations
	// **	
	double get_size() const ;
	void set_lai(double _l);
	void set_size(double _x, PlantTraits &traits);
	std::vector<double>::iterator set_state(std::vector<double>::iterator S, PlantTraits &traits);
	
	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - simulates growth over dt with constant assimilation rate A
	// ** 
	void grow_for_dt(double t, double dt, double &prod, double &litter_pool, double A, PlantTraits &traits);

};


} // namespace plant

#endif


