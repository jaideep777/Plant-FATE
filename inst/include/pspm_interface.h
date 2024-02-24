#ifndef PLANT_FATE_PSPM_INTERFACE_H_
#define PLANT_FATE_PSPM_INTERFACE_H_

#include <individual_base.h>
#include <solver.h>
#include "light_environment.h"
#include "climate.h"
#include "plant.h"

/// @defgroup libpspm_interface PSPM Interface
/// @brief    This is a collection of classes and functions used to interface with the PSPM Solver.

/// @defgroup ppa_module PPA
/// @brief    This module collects functions and classes that together implement the Perfect Plasticity Approximation 

#define STATE_DIM 1

/// @ingroup  libpspm_interface
/// @brief    This class entends the Plant class to interface with the PSPM Solver.
class PSPM_Plant : 
	public plant::Plant, public IndividualBase<STATE_DIM> 
{
	public:
	
	double t_birth = 0;
	
	/// @ingroup trait_evolution
	/// @brief   Variable names to print in the header corresponding to the output of PSPM_Plant::print.
	// FIXME: Make this static?
	std::vector<std::string> varnames = {"name", "|lma|", "| WD |", "D", "g", "lai", "mort", "seeds", "a", "c"}; // header corresponding to the print function below

	std::vector<std::string> statevarnames = {"lai", "mort"};  // header corresponding to state output (not used currently)

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	PSPM_Plant(); 

	void set_size(const std::array<double,STATE_DIM>& _x) override;

	double init_density(void * _env, double input_seed_rain) override;

	/// @ingroup libpspm_interface
	void preCompute(double t, void * _env) override;
	void afterStep(double t, void * _env);

	/// @addtogroup libpspm_interface
	/// Core demographic rate functions are specifically highlighted here
	/// even though they are members of the PSPM_Plant class.
	/// @{
	double establishmentProbability(double t, void * _env) override;

	std::array<double,STATE_DIM> growthRate(double t, void * _env) override;
	double mortalityRate(double t, void * _env) override;
	double birthRate(double t, void * _env) override;
	/// @}
	
	void init_accumulators(double t, void * _env) override;

	std::vector<double>::iterator set_accumulators(std::vector<double>::iterator &it) override;

	std::vector<double>::iterator get_accumulators(std::vector<double>::iterator &it) override;
	std::vector<double>::iterator get_accumulatorRates(std::vector<double>::iterator &it) override;

	void print(std::ostream &out = std::cout) const override;

	void save(std::ostream &fout) override;
	void restore(std::istream &fin) override;
};


/// @ingroup  libpspm_interface
/// @brief    This class provides a common interface to different components of the environment
class PSPM_Environment : 
	public env::LightEnvironment, public env::Climate
{
	public:
	void print(double t) override;
};


#endif
