#ifndef PLANT_FATE_PLANT_GEOMETRY_H_
#define PLANT_FATE_PLANT_GEOMETRY_H_

#include <cmath>

#include "plant_params.h"

namespace plant{

/// \ingroup physiology
class PlantGeometry{
	public:
	struct{
		// geometry traits
		double m, n;    ///< crown shape paramaters
		double a;       ///< height-diameter allometry
		double c;       ///< crown area allometry
		double fg;      ///< upper canopy gap fraction
		
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
	double lai;          ///< Crown leaf area index 
	double diameter;     ///< basal diameter (diameter at ground level)

	// variables calculated from state variables
	double height;                       ///< Plant height
	double crown_area;                   ///< Crown area
	double sapwood_fraction;             ///< Fraction of stem that is sapwood
	double functional_xylem_fraction;    ///< Fraction of funcitonal xylem in sapwood
	double rooting_depth;                ///< Rooting depth, calculated from coarse root biomass

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
	
	double diameter_at_height(double z, PlantTraits &traits);
	
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
	double leaf_mass(const PlantTraits &traits) const;
	double root_mass(const PlantTraits &traits) const;
	double sapwood_mass(const PlantTraits &traits) const;
	double sapwood_mass_real(const PlantTraits &traits) const;
	double stem_mass(const PlantTraits &traits) const;
	double coarse_root_mass(const PlantTraits &traits) const;
	double heartwood_mass(const PlantTraits &traits) const;
	double total_mass(const PlantTraits &traits) const;

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


