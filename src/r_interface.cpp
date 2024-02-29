#include <Rcpp.h>
using namespace Rcpp;

#include "traits_params.h"

RCPP_EXPOSED_CLASS(plant::PlantTraits)
RCPP_EXPOSED_CLASS(plant::PlantParameters)

#include "climate.h"
#include "light_environment.h"

RCPP_EXPOSED_CLASS(env::Clim);
RCPP_EXPOSED_CLASS(env::Climate);
RCPP_EXPOSED_CLASS(env::LightEnvironment);

#include "pspm_interface.h"

RCPP_EXPOSED_CLASS(EnvironmentBase);

#include "treelife.h"

RCPP_EXPOSED_CLASS(ErgodicEnvironment);

#include "plantfate.h"

RCPP_EXPOSED_CLASS(PSPM_Dynamic_Environment);


RCPP_MODULE(plantfate_module){
	class_ <plant::PlantTraits>("PlantTraits")
		.field("lma",          &plant::PlantTraits::lma)
		.field("zeta",         &plant::PlantTraits::zeta)
		.field("fcr",          &plant::PlantTraits::fcr)
		.field("hmat",         &plant::PlantTraits::hmat)
		.field("fhmat",        &plant::PlantTraits::fhmat)
		.field("seed_mass",    &plant::PlantTraits::seed_mass)
		.field("wood_density", &plant::PlantTraits::wood_density)
		.field("p50_xylem",    &plant::PlantTraits::p50_xylem)
		.field("K_leaf",       &plant::PlantTraits::K_leaf)
		.field("K_xylem",      &plant::PlantTraits::K_xylem)
		.field("b_leaf",       &plant::PlantTraits::b_leaf)
		.field("b_xylem",      &plant::PlantTraits::b_xylem)
		.field("m",            &plant::PlantTraits::m)
		.field("n",            &plant::PlantTraits::n)
		.field("a",            &plant::PlantTraits::a)
		.field("c",            &plant::PlantTraits::c)

		.constructor()
		.method("print", &plant::PlantTraits::print)
	;

	class_ <plant::PlantParameters>("PlantParameters")
		.field("cD0",          &plant::PlantParameters::cD0)
		.field("eD0",          &plant::PlantParameters::eD0)
		.field("cD1",          &plant::PlantParameters::cD1)
		.field("m_alpha",      &plant::PlantParameters::m_alpha)
		.field("m_beta",       &plant::PlantParameters::m_beta)
		.field("m_gamma",      &plant::PlantParameters::m_gamma)
		.field("eWD_alpha",    &plant::PlantParameters::eWD_alpha)
		.field("eWD_gamma",    &plant::PlantParameters::eWD_gamma)
		.field("cWD0",         &plant::PlantParameters::cWD0)
		.field("eWD",          &plant::PlantParameters::eWD)

		.constructor()
		.method("print", &plant::PlantParameters::print)	
	;

	class_ <env::LightEnvironment>("LightEnvironment")
		.constructor()
		.field("z_star", &env::LightEnvironment::z_star)
		.field("canopy_openness", &env::LightEnvironment::canopy_openness)
		.method("print", &env::LightEnvironment::print)
	;

	class_ <env::Clim>("Clim")
		.default_constructor()
		
		.field("tc", &env::Clim::tc)
		.field("ppfd_max", &env::Clim::ppfd_max)
		.field("ppfd", &env::Clim::ppfd)
		.field("vpd", &env::Clim::vpd)
		.field("co2", &env::Clim::co2)
		.field("elv", &env::Clim::elv)
		.field("swp", &env::Clim::swp)
	;

	class_ <env::Climate>("Climate")
		.default_constructor()

		.field("clim", &env::Climate::clim)

		// Dont expose these, because get/set functions in TreeLife and PlantFATE set both filename and update-flag
		// .field("metFile", &env::Climate::metFile)
		// .field("co2File", &env::Climate::co2File)
		// .field("update_met", &env::Climate::update_met)
		// .field("update_co2", &env::Climate::update_co2)

		.method("set", &env::Climate::set)
		.method("print", &env::Climate::print)
	;

	class_ <ErgodicEnvironment>("ErgodicEnvironment")
		.constructor()

		.derives<env::LightEnvironment>("LightEnvironment")
		.derives<env::Climate>("Climate")

		.method("print", &ErgodicEnvironment::print)
	;

	class_ <EnvironmentBase>("EnvironmentBase")
		// We do not want to export constructor to environment base because we dont want users to instantiate this
	;

	class_ <PSPM_Dynamic_Environment>("PSPM_Dynamic_Environment")
		// .derives<EnvironmentBase>("EnvironmentBase")
		.derives<env::LightEnvironment>("LightEnvironment")
		.derives<env::Climate>("Climate")

		.default_constructor()
	;


	class_ <LifeHistoryOptimizer>("LifeHistoryOptimizer")
		//.field("params_file", &LifeHistoryOptimizer::params_file)
		.field("env", &LifeHistoryOptimizer::C)
		.field("traits0", &LifeHistoryOptimizer::traits0)
		.field("par0", &LifeHistoryOptimizer::par0)
		
		.constructor<std::string>()
		.method("set_metFile", &LifeHistoryOptimizer::set_metFile)
		.method("set_co2File", &LifeHistoryOptimizer::set_co2File)

		.method("get_header", &LifeHistoryOptimizer::get_header)
		.method("get_state", &LifeHistoryOptimizer::get_state)

		.method("set_traits", &LifeHistoryOptimizer::set_traits)
		.method("get_traits", &LifeHistoryOptimizer::get_traits)
		.method("init", &LifeHistoryOptimizer::init)
		.method("printMeta", &LifeHistoryOptimizer::printMeta)
		.method("calcFitness", &LifeHistoryOptimizer::calcFitness)

		.method("grow_for_dt", &LifeHistoryOptimizer::grow_for_dt)
	;

	class_ <Patch>("Patch")
		.constructor<std::string>()
		.method("set_metFile", &Patch::set_metFile)
		.method("set_co2File", &Patch::set_co2File)
		.method("init", &Patch::init)
		.method("simulate", &Patch::simulate)
		.method("close", &Patch::close)

		.field("E", &Patch::E)

		.field("paramsFile", &Patch::paramsFile)
		.field("parent_dir", &Patch::parent_dir)
		.field("expt_dir",   &Patch::expt_dir)
		.field("save_state", &Patch::save_state)
		.field("state_outfile", &Patch::state_outfile)
		.field("config_outfile", &Patch::config_outfile)
		.field("continueFrom_stateFile", &Patch::continueFrom_stateFile)
		.field("continueFrom_configFile", &Patch::continueFrom_configFile)
		.field("continuePrevious", &Patch::continuePrevious)
		.field("evolve_traits", &Patch::evolve_traits)
		.field("y0", &Patch::y0)
		.field("yf", &Patch::yf)
		.field("ye", &Patch::ye)

		.field("traits0", &Patch::traits0)
		.field("par0", &Patch::par0)
	;
}

