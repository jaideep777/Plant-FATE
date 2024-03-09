#ifndef PLANT_FATE_PFATE_PLANTFATE_CONFIG_H_
#define PLANT_FATE_PFATE_PLANTFATE_CONFIG_H_

#include <vector>
#include <string>

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

	std::string time_unit;      ///< Time unit in which all time points and intervals are specified. E.g. days since 1850-01-01
	double timestep;            ///< Simulation timestep: Patches are stepped at this timestep and computations before/after timstep are done. Typically, these are climate update and seed-rain update 
	double T_cohort_insertion;  ///< Interval after which cohorts should be inserted, if using IEBT solver
	double T_seed_rain_avg;     ///< Interval over which seed rains should be averaged (multi-patch dynamics)
	double T_return;            ///< Return interval of disturbance (patch clearance)
	double T_invasion; 	        ///< Interval between successive species invasions

	double res; ///< initial resolution on size axis - remains constant for fixed-mesh methods

	std::string solver_method;

	std::vector<std::string> evolvable_traits;
	std::vector<double> trait_variances;
	std::vector<double> trait_scalars;
	double T_r0_avg;

};

} // namespace pfate

#endif