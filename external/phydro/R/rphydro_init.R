# rphydro_init.R

#' phydro
#'
#' phydro class
#' 
#' Imports
#' @useDynLib rphydro, .registration=T
#' @export phydro_module
#' @import Rcpp RcppEigen
"_PACKAGE"


Rcpp::loadModule(module="phydro_module", what=T) 


