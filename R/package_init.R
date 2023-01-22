# package_init.R

#' LifeHistoryOptimizer
#'
#' LifeHistoryOptimizer class
#'
#' Imports
#' @useDynLib PlantFATE, .registration=T
#' @export treelife_module, plantate_module
#' @import Rcpp
"_PACKAGE"


Rcpp::loadModule(module="treelife_module", what=T)
Rcpp::loadModule(module="plantfate_module", what=T)
