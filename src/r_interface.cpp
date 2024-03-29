#include <Rcpp.h>
using namespace Rcpp;

#include "treelife.h"

RCPP_MODULE(treelife_module){
	class_ <LifeHistoryOptimizer>("LifeHistoryOptimizer")
		.field("params_file", &LifeHistoryOptimizer::params_file)
		
		.constructor()
		.method("set_traits", &LifeHistoryOptimizer::set_traits)
		.method("get_traits", &LifeHistoryOptimizer::get_traits)
		.method("init", &LifeHistoryOptimizer::init)
		.method("printPlant", &LifeHistoryOptimizer::printPlant)
		.method("calcFitness", &LifeHistoryOptimizer::calcFitness)

		.method("grow_for_dt", &LifeHistoryOptimizer::grow_for_dt)
	;
}

#include "plantfate.h"

RCPP_MODULE(plantfate_module){
	class_ <Simulator>("Simulator")
		.constructor<std::string>()
		.method("init", &Simulator::init)
		.method("simulate", &Simulator::simulate)
		.method("close", &Simulator::close)

		.field("paramsFile", &Simulator::paramsFile)
		.field("parent_dir", &Simulator::parent_dir)
		.field("expt_dir", &Simulator::expt_dir)
		.field("met_file", &Simulator::met_file)
		.field("co2_file", &Simulator::co2_file)
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

