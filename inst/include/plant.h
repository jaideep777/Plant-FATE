#ifndef PLANT_FATE_PLANT_PLANT_H_
#define PLANT_FATE_PLANT_PLANT_H_
#include <fstream>
#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "utils/rk4.h"
#include "utils/moving_average.h"
#include "utils/initializer_v2.h"

namespace plant{

/// @defgroup physiology Physiology
/// @brief    This is a collection of classes which implement the core physiology of the plants in Plant-FATE.

/// @brief   The class for an individual plant
/// @ingroup physiology
class Plant{
	public:

	/// @brief Residual variables integrated via the ODE solver go in here. 
	struct{
		//double lai;         // these are in geometry 
		//double size;
		double mortality = 0;     ///< Cummulative mortality
		//double seed_pool = 0;
	} state;
	
	/// @brief Core demographic rates
	struct{
		double dlai_dt;    ///< Rate of change of LAI
		double dsize_dt;   ///< Growth rate
		double dmort_dt;   ///< Mortality rate
		double dseeds_dt;  ///< Fecundity rate 
		double rgr;        ///< Relative growth rate (RGR)
		// double dseeds_dt_pool;
		// double dseeds_dt_germ;
	} rates;	
		
	/// @brief Derivatives used for biomass partitioning
	struct {
		double dmass_dt_lai;
		double dmass_dt_rep;
		double dmass_dt_growth;
		double dmass_dt_lit;
		double dmass_dt_tot;
	} bp;

	/// @brief Result of whole-plant assimilation calculation, returned by Assimilator::calc_plant_assimilation_rate()
	PlantAssimilationResult res;

	public:
	//std::ofstream fmuh; // Cannot use streams here because we need copy-constructor for Plants, which in turn would need a copy constructor for streams, which is deleted.
	PlantTraits traits;   ///< Collection of all functional traits
	PlantParameters par;  ///< Collection of all model parameters that are not traits

	Assimilator assimilator; 
	PlantGeometry geometry;
	
	public:
		/// @brief  This function initializes the plant (traits, par, and geometry) from params and traits objects
		void init(const PlantParameters &_par, const PlantTraits &_traits);

		// /// @brief  This function initializes the plant (traits, par, and geometry) from an Initialzer object
		// void init(io::Initializer &I);

		// /// @brief  This function initializes the plant (traits, par, and geometry) from an ini file
		// void initFromFile(std::string file);

		/// @brief Set traits that are calculated from other traits (e.g., leaf_p50, a, c)
		void coordinateTraits();

		/// @addtogroup trait_evolution
		/// @{
		/// @brief Set values for evolvable traits from vector
		void set_evolvableTraits(std::vector<double> tvec);
		/// @brief Return values of evolvable traits in a vector
		std::vector<double> get_evolvableTraits();
		/// @}

		/// @brief  Set size (diameter and all associated variables) from x
		void set_size(double x); // FIXME: Is this function really needed?!
		/// @brief  Get plant biomass
		double get_biomass() const; // FIXME: Is this function really needed?!

		/// @brief LAI model
		template <class Env>
		double lai_model(PlantAssimilationResult &res, double _dmass_dt_tot, Env &env);

		/// @brief  Partition total biomass dm_dt_tot into various carbon pools
		/// @param dm_dt_tot  Total biomass to partition
		/// @param dm_dt_lai  Biomass that goes into LAI increment
		template <class Env>
		void partition_biomass(double dm_dt_tot, double dm_dt_lai, Env &env);

		// Core demographic rates
		/// @addtogroup libpspm_interface
		/// @{
		template <class Env>
		double size_growth_rate(double _dmass_dt_growth, Env &env);

		template <class Env>
		double mortality_rate(Env &env, double t);

		template <class Env>
		double fecundity_rate(double _dmass_dt_rep, Env &env);

		template <class Env>
		void calc_demographic_rates(Env &env, double t);
		/// @}

		/// @brief  Probability of survival during germination (i.e. until recruitment stage)
		template <class Env>
		double p_survival_germination(Env &env);

		/// @brief  Probability of survival during dispersal
		template <class Env>
		double p_survival_dispersal(Env &env);

		void print();

};


}	// namespace plant

#include "plant.tpp"

#endif
