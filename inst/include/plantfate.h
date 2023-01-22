#ifndef PLANT_FATE_PLANTFATE_H_
#define PLANT_FATE_PLANTFATE_H_

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
#include "trait_evolution.h"
#include "state_restore.h"

class Simulator{
	public:
	std::string paramsFile;
	std::string parent_dir, expt_dir, out_dir;

	std::string met_file;
	std::string co2_file;

	bool        save_state;
	std::string state_outfile;
	std::string config_outfile;

	std::string continueFrom_stateFile;
	std::string continueFrom_configFile;
	bool        continuePrevious;

	bool        evolve_traits;

	// Set up simulation start and end points
	double      y0;
	double      yf;
	double      ye;  // year in which trait evolution starts (need to allow this period because r0 is averaged over previous time)

	double t_clear = 105000;
	// t is years since 2000-01-01
	double delta_T;
	double timestep;

	std::string solver_method;

	io::Initializer          I;
	Solver                   S;
	PSPM_Dynamic_Environment E;

	SolverIO      sio;
	SpeciesProps  cwm;
	EmergentProps props; 

	public:
	Simulator(std::string params_file);
	
	void init(double tstart, double tend);

	void simulate();

	void close();

	private: 
	double runif(double rmin=0, double rmax=1);

	/// @brief     Calculate seed output of all species
	/// @param t   Current time 
	/// @param S   Solver
	/// @ingroup   trait_evolution
	/// @details   Species seed output rate is defined as,  
	///            \f[S = \int_{x_b}^{x_m}{f(s)u(s)ds}\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species) 
	void calc_seed_output(double t, Solver& S);


	// FIXME: Setting const input seed rain for mutants doesnt work. Is that a problem? 
	/// @brief     Calculate growth rates of all species and update seed input
	/// @param t   Current time
	/// @param dt  timestep over which growth rate is to be calculated
	/// @param S   Solver
	/// @ingroup   trait_evolution
	/// @details   Species growth rate is defined from the seed perspective, i.e.,
	///            \f[r = \frac{1}{\Delta t}log\left(\frac{S_\text{out}}{S_\text{in}}\right),\f] where \f$S\f$ is the seed rain (rate of seed production summed over all individuals of the species)
	void calc_r0(double t, double dt, Solver &S);

	void removeSpeciesAndProbes(Solver* S, MySpecies<PSPM_Plant>* spp);

	void addSpeciesAndProbes(Solver *S, std::string params_file, io::Initializer &I, double t, std::string species_name, double lma, double wood_density, double hmat, double p50_xylem);
	


};

#endif
