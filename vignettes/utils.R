#### Wrapper for classic pmodel with modified units ####

# P-model as per Wang et al 2017
# This needs PPFD in mol/m2/day

pmodel_wang17 = function(tc, ppfd, vpd, co2, elv, fapar, kphio, ...){
  
  out_analytical <- rpmodel::rpmodel(
    tc             = tc,
    vpd            = vpd,
    co2            = co2,
    elv            = elv,
    kphio          = kphio,
    beta           = 146,
    fapar          = fapar,
    ppfd           = ppfd*1e-6*86400, # p-model requires in mol/m2/day
    ...
  )
  
  # Convert some outputs to facilitate comparison
  return(list_modify(out_analytical,
                     gs = out_analytical$gs/86400*rpmodel::calc_patm(elv), # mol/m2/day/Pa --> mol/m2/s
                     gpp = out_analytical$gpp/1.03772448,  # gC/m2/day --> umol/m2/s
                     vcmax = out_analytical$vcmax/0.0864   # mol/m2/day --> umol/m2/s
  )
  )
}




#### Convenience Functions ####

get_std_photosythnesis_params = function(){ 
  ## Set P-model parameters
  kphio <- 0.05        # quantum yield efficiency
  # c_molmass <- 12.0107 # molar mass, g / mol
  
  ## Define environmental conditions
  tc <- 25             # temperature, deg C
  ppfd <- 400          # umol/m2/s
  # vpd  <- 1000         # Pa
  co2  <- 400          # ppm
  elv  <- 0            # m.a.s.l.
  fapar <- 0.7         # fraction
  
  p = rpmodel::calc_patm(elv)
  return (list(
    kmm = rpmodel::calc_kmm(tc, p),  # Why does this use std. atm pressure, and not p(z)?
    gammastar = rpmodel::calc_gammastar(tc, p),
    phi0 = kphio*rpmodel::calc_ftemp_kphio(tc),
    Iabs = ppfd*fapar,
    ca = co2*p*1e-6,  # Convert to partial pressure
    patm = p,
    delta = 0.00
  ))
}


# Combine Hydraulic and classic results

pmodel_calibrate_analytical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost = NULL, opt_hypothesis = "PM"){
  out_hydraulics = rphydro_analytical(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant, par_cost)
  
  out_analytical = pmodel_wang17(tc, ppfd, vpd, co2, elv, fapar, kphio)
  
  return(list(out_hydraulics = out_hydraulics,
              out_analytical = out_analytical))
}

pmodel_calibrate_numerical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost = NULL, opt_hypothesis = "PM"){
  out_hydraulics = rphydro_numerical(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant, par_cost)
  
  out_analytical = pmodel_wang17(tc, ppfd, vpd, co2, elv, fapar, kphio)
  
  return(list(out_hydraulics = out_hydraulics,
              out_analytical = out_analytical))
}

pmodel_calibrate_inst <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost = NULL, jmax, vcmax, opt_hypothesis = "PM"){
  if (is.null(par_cost)) par_cost = list(alpha=0.1, gamma=1)
  
  out_hydraulics = rphydro_instantaneous_analytical(vcmax, jmax, tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant, par_cost)
  
  out_analytical = pmodel_wang17(tc, ppfd, vpd, co2, elv, fapar, kphio)
  
  return(list(out_hydraulics = out_hydraulics,
              out_analytical = out_analytical))
}

