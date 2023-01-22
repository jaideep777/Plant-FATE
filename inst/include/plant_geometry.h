#ifndef PLANT_FATE_PLANT_GEOMETRY_H_
#define PLANT_FATE_PLANT_GEOMETRY_H_

#define _USE_MATH_DEFINES
#include <cmath>

#include "plant_params.h"

namespace plant{

/// @brief   Everything about the plant's physical dimensions
/// @ingroup physiology
/// @details Two variables completely define dimensional state of the 
///          plant: diameter and crown LAI. All other dimensions (height, 
///          crown area, leaf area, biomass, etc) are calculated from these two.
class PlantGeometry{
	public:
	/// @brief Traits that control dimensional scaling
	struct{
		// geometry traits
		double m, n;    ///< crown shape paramaters
		double a;       ///< height-diameter allometry
		double c;       ///< crown area allometry
		double fg;      ///< upper canopy gap fraction
		
		// Precomputed Geometric parameters
		double eta_c;   ///< Stem taper coefficient
		double pic_4a;  ///< \f$\pi c / (4a)\f$
		double zm_H;    ///< \f$z_m / H\f$
		double qm;      ///< \f$q(z_m)\f$

		// allocation
		double dmat;    ///< Diameter at reproductive maturity, calculated as diameter when H = `fhmat x` hmat
	} geom;

	public:
	// current state
	double lai;          ///< Crown leaf area index 
	double diameter;     ///< basal diameter (diameter at ground level)

	// variables calculated from state variables
	double height;                       ///< Plant height
	double crown_area;                   ///< Crown area
	double sapwood_fraction;             ///< Fraction of stem cross sectional area that is sapwood
	double functional_xylem_fraction;    ///< Fraction of funcitonal xylem in sapwood
	double rooting_depth;                ///< Rooting depth, calculated from coarse root biomass

	// ode-based calculations of sapwood and heartwood (for debug)
	double sap_frac_ode = 1;
	double sapwood_mass_ode = 0;
	double heart_mass_ode = 0;
	double k_sap;

	public:

	/// @brief  Initialize geometry from traits, precompute any necessary variables
	void init(PlantParameters &par, PlantTraits &traits);


	/// @brief  The height at which crown radius is maximum.
	double zm();


	/// @brief Vertical profiles of crown and stem  
	/// @{
	double q(double z);
	double diameter_at_height(double z, PlantTraits &traits);
	/// @brief  Potential crown projection area at height z. 
	double crown_area_extent_projected(double z, PlantTraits &traits);
	/// @brief  Realized crown projection area at height z. 
	double crown_area_above(double z, PlantTraits &traits);
	/// @}


	/// @brief Derivatives required for biomass partitioning
	/// @{  
	double dsize_dmass(PlantTraits &traits) const ;
	double dreproduction_dmass(PlantParameters &par, PlantTraits &traits);
	/// @}


	/// Rate of change of leaf mass due to change in LAI
	double dmass_dt_lai(double &dL_dt, double dmass_dt_max, PlantTraits &traits);


	/// @brief Get biomass in various carbon pools.
	/// @{
	double leaf_mass(const PlantTraits &traits) const;
	double root_mass(const PlantTraits &traits) const;
	double sapwood_mass(const PlantTraits &traits) const;
	double sapwood_mass_real(const PlantTraits &traits) const;
	double stem_mass(const PlantTraits &traits) const;
	double coarse_root_mass(const PlantTraits &traits) const;
	double heartwood_mass(const PlantTraits &traits) const;
	double total_mass(const PlantTraits &traits) const;
	/// @}


	// These functions are used to get and set state variables
	/// @{
	/// Get size (diameter) - does not alter state
	double get_size() const ;
	/// Set the crown LAI and properties that change with LAI
	void set_lai(double _l);
	/// Set plant size (diameter) and other variables that scale with size  
	void set_size(double _x, PlantTraits &traits);
	/// Set size and lai, the two state variables that define plant geometry
	std::vector<double>::iterator set_state(std::vector<double>::iterator S, PlantTraits &traits);
	/// @}


	// ** 
	// ** Simple growth simulator for testing purposes
	// ** - simulates growth over dt with constant assimilation rate A
	// ** 
	void grow_for_dt(double t, double dt, double &prod, double &litter_pool, double A, PlantTraits &traits);

};


} // namespace plant

#endif


