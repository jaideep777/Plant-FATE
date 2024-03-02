// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

//#include <memory> // std::shared_ptr
//#include <vector>
//#include <plant/ode_interface.h>
//#include <plant/plant_internals.h>
#include <limits>
#include <iostream>
#include <vector>

//#include "environment.h"
#include "pn_integrator.h"

namespace plant {


struct PlantParameters{
	// * Individual allometry
	// Canopy shape parameters
	double eta, eta_c;
	// Sapwood area per leaf area
	double theta;
	// Empirical constants for scaling relationships
	double a_l1, a_l2, a_r1;
	// Bark area per sapwood area
	double a_b1;
	// * Production
	// Respiration constants
	double r_s, r_b, r_r, r_l;
	// Yield = carbon fixed in tissue per carbon assimilated;
	double a_y;
	// Conversion factor
	double a_bio;
	// Leaf, bark sapwood, and root turnover rates
	double k_l, k_b, k_s, k_r;
	// Leaf productivity parameters  - only used when no N reallocation
	double a_p1, a_p2;

	// * Seed production
	// Accessory cost of reproduction, per seed
	double a_f3;
	// Proportion production allocated to reproduction
	double a_f1;
	// Size range across which individuals mature
	double a_f2;

	// * Mortality
	// Probability of survival during dispersal
	double S_D;
	// Parameter for seedling mortality
	double a_d0;
	// Baseline structural mortality rate
	double d_I;
	double c_d0;
	double c_d1;
	// Baseline for growth mortality rate
	double a_dG1;
	// Coefficient for dry mass production in mortality function
	double a_dG2;

	//  // Height and leaf area of a (germinated) seed
	double height_0;
	double area_leaf_0;
};
//void initPlantParameters(PlantParameters &par);



//extern PlantParameters par;	// need this as global variable because in GPU implementation this will be __constant__

class Plant {
	public:
	
	// * Core traits
	double lma       = 0.0825; //0.1978791;  // Leaf mass per area [kg / m2]
	double rho       = 608.0;      // Wood density [kg/m3]
	double hmat      = 16.5958691; // Height at maturation [m]
	double omega     = 3.8e-5;     // Seed mass [kg]

	double trait_eta = 12.0;

	PlantParameters par;	// JAI: making this class member again for multispecies version

	struct{
		double height            = std::numeric_limits<double>::quiet_NaN(); 
		double height_dt         = std::numeric_limits<double>::quiet_NaN();
		double area_leaf         = 0.0;
		double mortality         = 0.0;
		double mortality_dt      = std::numeric_limits<double>::quiet_NaN();
		double fecundity         = 0.0;
		double fecundity_dt      = std::numeric_limits<double>::quiet_NaN();
		double area_heartwood    = 0.0;
		double area_heartwood_dt = std::numeric_limits<double>::quiet_NaN();
		double mass_heartwood    = 0.0;
		double mass_heartwood_dt = std::numeric_limits<double>::quiet_NaN();
	} vars;	// TODO: Remove this eventually 
	
	static pn::Integrator plantIntegrator;

	public:
	
	Plant();

	void initParameters();

	friend std::ostream& operator<<(std::ostream& os, const Plant& p);

	
	double height() const;
	double height_dt() const;
	void   set_height(double x);

	double mortality() const;
	double mortality_dt() const;
	void   set_mortality(double x);

	double fecundity() const;
	double fecundity_dt() const;
	void   set_fecundity(double x);

	double area_heartwood() const;
	double area_heartwood_dt() const;
	void   set_area_heartwood(double x);

	double mass_heartwood() const;
	double mass_heartwood_dt() const;
	void   set_mass_heartwood(double x);


  // [eqn 2] area_leaf (inverse of [eqn 3])
  double area_leaf(double height) const;

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  double mass_leaf(double area_leaf) const;

  // [eqn 4] area and mass of sapwood
  double area_sapwood(double area_leaf) const;
  double mass_sapwood(double area_sapwood, double height) const;

  // [eqn 5] area and mass of bark
  double area_bark(double area_leaf) const;
  double mass_bark (double area_bark, double height) const;

  double area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const;
  double diameter_stem(double area_stem) const;

  // [eqn 7] Mass of (fine) roots
  double mass_root(double area_leaf) const;

  // [eqn 8] Total Mass
  double mass_live(double mass_leaf, double mass_bark,
                   double mass_sapwood, double mass_root) const;

  double mass_total(double mass_leaf, double mass_bark, double mass_sapwood,
                    double mass_heartwood, double mass_root) const;

  double mass_above_ground(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const;

  template <class Environment>
	  void compute_vars_phys(const Environment& environment, bool reuse_intervals=false);

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  template <class Environment>
	  double assimilation(const Environment& environment, double height,
                      double area_leaf, bool reuse_intervals);
  // Used internally, corresponding to the inner term in [eqn 12]
  //double compute_assimilation_x(double x, double height,
  //                              const Environment& environment) const;
  //double compute_assimilation_h(double h, double height,
  //                              const Environment& environment) const;
  //double compute_assimilation_p(double p, double height,
  //                              const Environment& environment) const;
  // [Appendix S6] Per-leaf photosynthetic rate.
  double assimilation_leaf(double x) const;



  // [eqn 13] Total maintenance respiration
  double respiration(double mass_leaf, double mass_sapwood,
                     double mass_bark, double mass_root) const;

  double respiration_leaf(double mass) const;
  double respiration_bark(double mass) const;
  double respiration_sapwood(double mass) const;
  double respiration_root(double mass) const;

  // [eqn 14] Total turnover
  double turnover(double mass_leaf, double mass_bark,
                  double mass_sapwood, double mass_root) const;
  double turnover_leaf(double mass) const;
  double turnover_bark(double mass) const;
  double turnover_sapwood(double mass) const;
  double turnover_root(double mass) const;

  // [eqn 15] Net production
  //double net_mass_production_dt_A(double assimilation, double respiration,
  //                                double turnover) const;
  template <class Environment>
	  double net_mass_production_dt(const Environment& environment,
                                double height, double area_leaf_,
                                bool reuse_intervals=false);

  // [eqn 16] Fraction of whole plan growth that is leaf
  double fraction_allocation_reproduction(double height) const;
  double fraction_allocation_growth(double height) const;
  // [eqn 17] Rate of offspring production
  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  double darea_leaf_dmass_live(double area_leaf) const;

  // change in height per change in leaf area
  double dheight_darea_leaf(double area_leaf) const;
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  double dmass_leaf_darea_leaf(double area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  double dmass_sapwood_darea_leaf(double area_leaf) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  double dmass_bark_darea_leaf(double area_leaf) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  double dmass_root_darea_leaf(double area_leaf) const;
  // Growth rate of basal diameter_stem per unit stem area
  double ddiameter_stem_darea_stem(double area_stem) const;
  // Growth rate of components per unit time:
  double area_leaf_dt(double area_leaf_dt) const;
  double area_sapwood_dt(double area_leaf_dt) const;
  double area_heartwood_dt(double area_leaf) const;
  double area_bark_dt(double area_leaf_dt) const;
  double area_stem_dt(double area_leaf, double area_leaf_dt) const;
  double diameter_stem_dt(double area_stem, double area_stem_dt) const;
  double mass_root_dt(double area_leaf,
                       double area_leaf_dt) const;
  double mass_live_dt(double fraction_allocation_reproduction,
                       double net_mass_production_dt) const;
  double mass_total_dt(double fraction_allocation_reproduction,
                        double net_mass_production_dt,
                        double mass_heartwood_dt) const;
  double mass_above_ground_dt(double area_leaf,
                               double fraction_allocation_reproduction,
                               double net_mass_production_dt,
                               double mass_heartwood_dt,
                               double area_leaf_dt) const;

  double mass_heartwood_dt(double mass_sapwood) const;

  double mass_live_given_height(double height) const;
  double height_given_mass_leaf(double mass_leaf_) const;


  double mortality_dt(double productivity_area, double cumulative_mortality) const;
  double mortality_growth_independent_dt()const ;
  double mortality_growth_dependent_dt(double productivity_area) const;
  // [eqn 20] Survival of seedlings during germination
  template <class Environment>
	  double germination_probability(const Environment& environment);

  // * Competitive environment
  // [eqn 11] total leaf area above height above height `z` for given plant
  double area_leaf_above(double z, double height, double area_leaf) const;
  // [eqn  9] Probability density of leaf area at height `z`
  double q(double z, double height) const;
  // [eqn 10] Fraction of leaf area above height `z`
  double Q(double z, double height) const;
  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  double Qp(double x, double height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;


	// New ODE interface
	void set_state(const double * full_state_array);
	std::vector <double> get_state();
	void get_rates(double * all_rates_array);


};

#include "../src/plant.tpp"

}

#endif
