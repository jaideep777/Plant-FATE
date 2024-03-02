# These are R wrappers for exposed Rcpp functions. These are 
# defined so that the functions can be used with named arguments

rr_calc_ftemp_inst_vcmax = function(tc, tg, tref, method){
  r_calc_ftemp_inst_vcmax(tc, tg, tref, method)
}

rr_calc_ftemp_inst_jmax = function(tc, tg, th, tref, method){
  r_calc_ftemp_inst_jmax(tc, tg, th, tref, method)
}

