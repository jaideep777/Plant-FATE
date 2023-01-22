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
	;
}

