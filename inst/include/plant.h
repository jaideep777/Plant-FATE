#ifndef PLANT_FATE_PLANT_PLANT_H_
#define PLANT_FATE_PLANT_PLANT_H_
#include <fstream>
#include "plant_params.h"
#include "plant_geometry.h"
#include "assimilation.h"
#include "utils/rk4.h"
#include "utils/moving_average.h"

namespace plant{

/// \defgroup physiology Physiology
/// \brief    This is a collection of classes which implement the core physiology of the plants in Plant-FATE.

/// \ingroup physiology
class Plant{
	public:
	// ** core state variables **
	// Integrated variables go in this struct, if not part of Plant already
	struct{
		//double lai;         // these are in geometry 
		//double size;
		double mortality = 0;     // cummulative mortality
		//double seed_pool = 0;
	} state;
	
	// ** core rates **
	struct{
		double dlai_dt;
		double dsize_dt;   // growth rate
		double dmort_dt;   // mortality rate
		double dseeds_dt;  // fecundity rate 
		// double dseeds_dt_pool;
		// double dseeds_dt_germ;
		double rgr;        // relative growth rate
	} rates;	
		
	// results of biomass partitioning
	struct {
		double dmass_dt_lai;
		double dmass_dt_rep;
		double dmass_dt_growth;
		double dmass_dt_lit;
		double dmass_dt_tot;
	} bp;


	PlantAssimilationResult res;

	// seed output history
	//MovingAverager seeds_hist;

	public:
	//std::ofstream fmuh; // Cannot use streams here because we need copy-constructor for Plants, which in turn would need a copy constructor for streams, which is deleted.
	PlantTraits traits;
	PlantParameters par;

	Assimilator assimilator; // to use pointers here, need to apply rule of 5
	PlantGeometry geometry;
	
	public:

	void initParamsFromFile(std::string file);
	void coordinateTraits();
	
	/// @addtogroup trait_evolution
	/// @{
	void set_evolvableTraits(std::vector<double> tvec);
	std::vector<double> get_evolvableTraits();
	/// @}

	void set_size(double x);

	double get_biomass() const;
	
	// LAI model
	template<class Env>
	double lai_model(PlantAssimilationResult& res, double _dmass_dt_tot, Env &env);

	template<class Env>
	void partition_biomass(double dm_dt_tot, double dm_dt_lai, Env &env);

	/// @addtogroup libpspm_interface
	/// @{
	// demographics
	template<class Env>
	double size_growth_rate(double _dmass_dt_growth, Env &env);

	template<class Env>
	double mortality_rate(Env &env, double t);

	template<class Env>
	double fecundity_rate(double _dmass_dt_rep, Env &env);

	template<class Env>
	void calc_demographic_rates(Env &env, double t);
	/// @}

	template<class Env>
	double p_survival_germination(Env &env);

	template<class Env>
	double p_survival_dispersal(Env &env);


	void print();

};


}	// namespace plant

#include "plant.tpp"

#endif
