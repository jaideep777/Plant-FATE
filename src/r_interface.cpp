#include <Rcpp.h>
using namespace Rcpp;

#include "climate.h"
#include "light_environment.h"
#include "treelife.h"

RCPP_EXPOSED_CLASS(ErgodicEnvironment);
RCPP_EXPOSED_CLASS(env::LightEnvironment);
RCPP_EXPOSED_CLASS(env::Climate);
RCPP_EXPOSED_CLASS(env::Clim);


RCPP_MODULE(treelife_module){
	class_ <env::LightEnvironment>("LightEnvironment")
		.constructor()
		.field("z_star", &env::LightEnvironment::z_star)
		.field("canopy_openness", &env::LightEnvironment::canopy_openness)
		.method("print", &env::LightEnvironment::print)
	;

	class_ <env::Clim>("Clim")
		.field("tc", &env::Clim::tc)
		.field("ppfd_max", &env::Clim::ppfd_max)
		.field("ppfd", &env::Clim::ppfd)
		.field("vpd", &env::Clim::vpd)
		.field("co2", &env::Clim::co2)
		.field("elv", &env::Clim::elv)
		.field("swp", &env::Clim::swp)
	;

	class_ <env::Climate>("Climate")
		.constructor()
		.field("clim", &env::Climate::clim)

		.field("metFile", &env::Climate::metFile)
		.field("co2File", &env::Climate::co2File)
		.field("update_met", &env::Climate::update_met)
		.field("update_co2", &env::Climate::update_co2)

		.method("set", &env::Climate::set)
		.method("print", &env::Climate::print)
	;

	class_ <ErgodicEnvironment>("ErgodicEnvironment")
		.constructor()
		.derives<env::LightEnvironment>("LightEnvironment")
		.derives<env::Climate>("Climate")

		.method("print", &ErgodicEnvironment::print)
	;

	class_ <LifeHistoryOptimizer>("LifeHistoryOptimizer")
		.field("params_file", &LifeHistoryOptimizer::params_file)
		.field("env", &LifeHistoryOptimizer::C)
		
		.constructor()
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
}

#include "plantfate.h"
#include "pspm_interface.h"

RCPP_EXPOSED_CLASS(EnvironmentBase);
RCPP_EXPOSED_CLASS(PSPM_Dynamic_Environment);

RCPP_MODULE(plantfate_module){
	class_ <EnvironmentBase>("EnvironmentBase")
	;

	class_ <PSPM_Dynamic_Environment>("PSPM_Dynamic_Environment")
	;

	class_ <Simulator>("Simulator")
		.constructor<std::string>()
		.method("set_metFile", &Simulator::set_metFile)
		.method("set_co2File", &Simulator::set_co2File)
		.method("init", &Simulator::init)
		.method("simulate", &Simulator::simulate)
		.method("close", &Simulator::close)

		.field("E", &Simulator::E)

		.field("paramsFile", &Simulator::paramsFile)
		.field("parent_dir", &Simulator::parent_dir)
		.field("expt_dir", &Simulator::expt_dir)
		.field("save_state", &Simulator::save_state)
		.field("state_outfile", &Simulator::state_outfile)
		.field("config_outfile", &Simulator::config_outfile)
		.field("continueFrom_stateFile", &Simulator::continueFrom_stateFile)
		.field("continueFrom_configFile", &Simulator::continueFrom_configFile)
		.field("continuePrevious", &Simulator::continuePrevious)
		.field("evolve_traits", &Simulator::evolve_traits)
		.field("y0", &Simulator::y0)
		.field("yf", &Simulator::yf)
		.field("ye", &Simulator::ye)
	;
}

