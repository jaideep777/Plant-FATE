#ifndef PLANT_FATE_PFATE_PATCH_H_
#define PLANT_FATE_PFATE_PATCH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include <solver.h>
#include "pspm_interface.h"
#include "trait_reader.h"
#include "community_properties.h"
#include "adaptive_species.h"
#include "state_restore.h"
#include "climate_stream.h"

namespace pfate{

class PlantFateConfig{

	public:
	std::string paramsFile;

	std::string out_dir;
	std::string parent_dir, expt_dir;

	// std::string i_met_file;
	// std::string a_met_file;
	// std::string co2_file;

	bool        save_state;
	std::string state_outfile;
	std::string config_outfile;

	std::string continueFrom_stateFile;
	std::string continueFrom_configFile;
	bool        continuePrevious;
	int         saveStateInterval;

	int         n_species;
	std::string traits_file;
	bool        evolve_traits;

	// Set up simulation start and end points
	double      y0;
	double      yf;
	double      ye;  // year in which trait evolution starts (need to allow this period because r0 is averaged over previous time)

	// t is years since 2000-01-01
	double timestep;            ///< Simulation timestep: Patches are stepped at this timestep and computations before/after timstep are done. Typically, these are climate update and seed-rain update 
	double T_cohort_insertion;  ///< Interval after which cohorts should be inserted, if using IEBT solver
	double T_seed_rain_avg;     ///< Interval over which seed rains should be averaged (multi-patch dynamics)
	double T_return;            ///< Return interval of disturbance (patch clearance)
	double T_invasion; 	        ///< Interval between successive species invasions

	double res; // initial resolution on size axis - remains constant for fixed-mesh methods

	std::string solver_method;

	std::vector<double> trait_variances;
	std::vector<double> trait_scalars;
	double T_r0_avg;

};


class Patch{
	public:
	PlantFateConfig config;

	double t_next_disturbance;
	double t_next_invasion;
	double t_next_savestate;
	double t_next_writestate;

	plant::PlantParameters par0;
	plant::PlantTraits traits0;

	env::ClimateStream  climate_stream;  // should be moved out of patch

	// io::Initializer     I;
	Solver              S;
	PSPM_Environment    E;

	SolverIO      sio;
	SpeciesProps  cwm;
	EmergentProps props; 

	public:
	Patch(std::string params_file);
	
	void set_i_metFile(std::string file);
	void set_a_metFile(std::string file);
	void set_co2File(std::string co2file);

	void init(double tstart, double tend);

	void simulate_to(double t);

	void update_climate(double julian_time, env::ClimateStream& c_stream);

	void update_climate(double co2, double tc, double vpd, double ppfd, double swp);
	void update_climate_acclim(double t_julian, double co2, double tc, double vpd, double ppfd, double swp);

	void simulate();

	void close();

	private: 
	double runif(double rmin=0, double rmax=1);

	// FIXME: Setting const input seed rain for mutants doesnt work. Is that a problem? 
	/// @brief     Calculate growth rates of all species and update seed input
	/// @param t   Current time
	/// @param dt  timestep over which growth rate is to be calculated
	/// @param S   Solver
	/// @ingroup   trait_evolution
	/// @details   Species growth rate is defined from the seed perspective, i.e.,
	///            \f[r = \frac{1}{\Delta t}log\left(\frac{S_\text{out}}{S_\text{in}}\right),\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species)
	void calc_seedrain_r0(double t);

	void removeSpeciesAndProbes(AdaptiveSpecies<PSPM_Plant>* spp);
	void addSpeciesAndProbes(double t, const plant::PlantTraits& traits);
	void shuffleSpecies();

	void removeDeadSpecies(double t);
	void addRandomSpecies(double t);
	void evolveTraits(double t, double dt_evolution);
	void disturbPatch(double t);
};

} // namespace pfate

#endif

