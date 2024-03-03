#include <Rcpp.h>
using namespace Rcpp;

#include "traits_params.h"

RCPP_EXPOSED_CLASS_NODECL(plant::PlantTraits)
RCPP_EXPOSED_CLASS_NODECL(plant::PlantParameters)

#include "climate.h"
#include "light_environment.h"

RCPP_EXPOSED_CLASS_NODECL(pfate::env::Clim);
RCPP_EXPOSED_CLASS_NODECL(pfate::env::Climate);
RCPP_EXPOSED_CLASS_NODECL(pfate::env::LightEnvironment);

#include "pspm_interface.h"

RCPP_EXPOSED_CLASS_NODECL(EnvironmentBase);
RCPP_EXPOSED_CLASS_NODECL(pfate::PSPM_Environment);

#include "life_history.h"

RCPP_EXPOSED_CLASS_NODECL(pfate::ErgodicEnvironment);

#include "plantfate_patch.h"

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
		.field("sm_xylem",      &plant::PlantTraits::sm_xylem)
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
		.field("m_hydraulic",  &plant::PlantParameters::m_hydraulic)

		.constructor()
		.method("print", &plant::PlantParameters::print)	
	;

	class_ <pfate::env::LightEnvironment>("LightEnvironment")
		.constructor()
		.field("z_star", &pfate::env::LightEnvironment::z_star)
		.field("canopy_openness", &pfate::env::LightEnvironment::canopy_openness)
		.method("print", &pfate::env::LightEnvironment::print)
	;

	class_ <pfate::env::Clim>("Clim")
		.default_constructor()
		
		.field("tc",    &pfate::env::Clim::tc)
		.field("ppfd",  &pfate::env::Clim::ppfd)
		.field("vpd",   &pfate::env::Clim::vpd)
		.field("co2",   &pfate::env::Clim::co2)
		.field("elv",   &pfate::env::Clim::elv)
		.field("swp",   &pfate::env::Clim::swp)
		.field("rn",    &pfate::env::Clim::rn)
		.field("vwind", &pfate::env::Clim::vwind)
	;

	class_ <pfate::env::Climate>("Climate")
		.default_constructor()

		.field("clim_inst", &pfate::env::Climate::clim_inst)
		.field("clim_acclim", &pfate::env::Climate::clim_acclim)

		// Dont expose these, because get/set functions in TreeLife and PlantFATE set both filename and update-flag
		// .field("metFile", &pfate::env::Climate::metFile)
		// .field("co2File", &pfate::env::Climate::co2File)
		// .field("update_met", &pfate::env::Climate::update_met)
		// .field("update_co2", &pfate::env::Climate::update_co2)

		// FIXME: These setters DO NOT WORK on the env object within LifeHstoryOptimizer, e.g., lho$env$init_co2()
		//        Hence commenting out. Use the overloads in LHO instead.
		// .method("init_co2", &pfate::env::Climate::init_co2)
		// .method("set_forcing_acclim", &pfate::env::Climate::set_forcing_acclim)

		.method("print", &pfate::env::Climate::print)
	;

	class_ <pfate::ErgodicEnvironment>("ErgodicEnvironment")
		.constructor()

		.derives<pfate::env::LightEnvironment>("LightEnvironment")
		.derives<pfate::env::Climate>("Climate")

		.method("print", &pfate::ErgodicEnvironment::print)
	;

	class_ <EnvironmentBase>("EnvironmentBase")
		// We do not want to export constructor to environment base because we dont want users to instantiate this
	;

	class_ <pfate::PSPM_Environment>("PSPM_Environment")
		// .derives<EnvironmentBase>("EnvironmentBase")
		.derives<pfate::env::LightEnvironment>("LightEnvironment")
		.derives<pfate::env::Climate>("Climate")

		.default_constructor()
	;


	class_ <pfate::LifeHistoryOptimizer>("LifeHistoryOptimizer")
		//.field("params_file", &pfate::LifeHistoryOptimizer::params_file)
		.field("env", &pfate::LifeHistoryOptimizer::C)
		.field("traits0", &pfate::LifeHistoryOptimizer::traits0)
		.field("par0", &pfate::LifeHistoryOptimizer::par0)
		
		.constructor<std::string>()
		.method("set_i_metFile", &pfate::LifeHistoryOptimizer::set_i_metFile)
		.method("set_a_metFile", &pfate::LifeHistoryOptimizer::set_a_metFile)
		.method("set_co2File", &pfate::LifeHistoryOptimizer::set_co2File)
		.method("init_co2", &pfate::LifeHistoryOptimizer::init_co2)

		.method("get_header", &pfate::LifeHistoryOptimizer::get_header)
		.method("get_state", &pfate::LifeHistoryOptimizer::get_state)

		// .method("set_traits", &pfate::LifeHistoryOptimizer::set_traits)
		// .method("get_traits", &pfate::LifeHistoryOptimizer::get_traits)
		.method("init", &pfate::LifeHistoryOptimizer::init)
		.method("printMeta", &pfate::LifeHistoryOptimizer::printMeta)
		.method("calcFitness", &pfate::LifeHistoryOptimizer::calcFitness)

		.method("grow_for_dt", &pfate::LifeHistoryOptimizer::grow_for_dt)
	;

	class_ <pfate::Patch>("Patch")
		.constructor<std::string>()
		.method("set_i_metFile", &pfate::Patch::set_i_metFile)
		.method("set_a_metFile", &pfate::Patch::set_a_metFile)
		.method("set_co2File", &pfate::Patch::set_co2File)
		.method("init", &pfate::Patch::init)
		.method("simulate", &pfate::Patch::simulate)
		.method("close", &pfate::Patch::close)

		.field("E", &pfate::Patch::E)

		.field("traits0", &pfate::Patch::traits0)
		.field("par0", &pfate::Patch::par0)
	;

	class_ <pfate::PlantFateConfig>("pfate::PlantFateConfig")
		.field("paramsFile", &pfate::PlantFateConfig::paramsFile)
		.field("parent_dir", &pfate::PlantFateConfig::parent_dir)
		.field("expt_dir",   &pfate::PlantFateConfig::expt_dir)
		.field("save_state", &pfate::PlantFateConfig::save_state)
		.field("state_outfile", &pfate::PlantFateConfig::state_outfile)
		.field("config_outfile", &pfate::PlantFateConfig::config_outfile)
		.field("continueFrom_stateFile", &pfate::PlantFateConfig::continueFrom_stateFile)
		.field("continueFrom_configFile", &pfate::PlantFateConfig::continueFrom_configFile)
		.field("continuePrevious", &pfate::PlantFateConfig::continuePrevious)
		.field("evolve_traits", &pfate::PlantFateConfig::evolve_traits)
		.field("y0", &pfate::PlantFateConfig::y0)
		.field("yf", &pfate::PlantFateConfig::yf)
		.field("ye", &pfate::PlantFateConfig::ye)
	;

}

