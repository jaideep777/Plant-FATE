#ifndef PLANT_FATE_PSPM_INTERFACE_H_
#define PLANT_FATE_PSPM_INTERFACE_H_

#include <solver.h>
#include "light_environment.h"
#include "climate.h"
#include "plant.h"

/// @defgroup libpspm_interface PSPM Interface
/// @brief    This is a collection of classes and functions used to interface with the PSPM Solver.

/// @defgroup ppa_module PPA
/// @brief    This module collects functions and classes that together implement the Perfect Plasticity Approximation 

/// @ingroup  libpspm_interface
/// @brief    This class entends the Plant class to interface with the PSPM Solver.
class PSPM_Plant : public plant::Plant {
	public:
	
	double t_birth = 0;
	
	/// @ingroup trait_evolution
	/// @brief   Variable names to print in the header corresponding to the output of PSPM_Plant::print.
	std::vector<std::string> varnames = {"name", "|lma|", "| WD |", "D", "g", "lai", "mort", "seeds"}; // header corresponding to the print function below

	std::vector<std::string> statevarnames = {"lai", "mort"};  // header corresponding to state output (not used currently)

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	PSPM_Plant(); 

	void set_size(double _x);

	double init_density(double x, void * _env, double input_seed_rain);

	/// @ingroup libpspm_interface
	void preCompute(double x, double t, void * _env);
	void afterStep(double x, double t, void * _env);

	/// @addtogroup libpspm_interface
	/// Core demographic rate functions are specifically highlighted here
	/// even though they are members of the PSPM_Plant class.
	/// @{
	double establishmentProbability(double t, void * _env);

	double growthRate(double x, double t, void * _env);
	double mortalityRate(double x, double t, void * _env);
	double birthRate(double x, double t, void * _env);
	/// @}
	
	void init_state(double t, void * _env);

	std::vector<double>::iterator set_state(std::vector<double>::iterator &it);

	std::vector<double>::iterator get_state(std::vector<double>::iterator &it);
	std::vector<double>::iterator get_rates(std::vector<double>::iterator &it);

	void print(std::ostream &out = std::cout);

	void save(std::ostream &fout);
	void restore(std::istream &fin);
};

/// @ingroup  libpspm_interface
/// @brief    Environment class for interfacing with the PSPM Solver
class PSPM_Dynamic_Environment : public EnvironmentBase, public env::LightEnvironment, public env::Climate{
	public:
	double projected_crown_area_above_z(double t, double z, Solver *S);
	double fapar_layer(double t, int layer, Solver *S);
	void computeEnv(double t, Solver *S, std::vector<double>::iterator _S, std::vector<double>::iterator _dSdt);
	void print(double t);
};





#endif
