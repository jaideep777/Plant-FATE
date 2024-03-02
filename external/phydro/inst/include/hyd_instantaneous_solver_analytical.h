#ifndef PHYDRO_INST_ANALYTICAL_SOLVER_H
#define PHYDRO_INST_ANALYTICAL_SOLVER_H

#include "hyd_transpiration.h"
#include "hyd_photosynthesis.h"

namespace phydro{

inline double calc_dP_ddpsi(double dpsi, double vcmax, double jmax, double psi_soil, ParPlant par_plant, ParEnv par_env, ParPhotosynth par_photosynth, ParCost par_cost){
	double gstar = par_photosynth.gammastar/par_photosynth.patm*1e6;
	double Km = par_photosynth.kmm/par_photosynth.patm*1e6;
	double ca = par_photosynth.ca/par_photosynth.patm*1e6;
	double br = par_photosynth.delta;
	double y = par_cost.gamma;

	double Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env);
	double gs = calc_gs_from_Q(Q, psi_soil, par_plant, par_env);
	auto Assim = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth);
	double P = Assim.a - y*dpsi*dpsi;

	double dpsi1 = dpsi+1e-6;
	double Q1 = calc_sapflux(dpsi1, psi_soil, par_plant, par_env);
	double gs1 = calc_gs_from_Q(Q1, psi_soil, par_plant, par_env);
	auto Assim1 = calc_assimilation_limiting(vcmax, jmax, gs1, par_photosynth);
	double P1 = Assim1.a - y*(dpsi1)*(dpsi1);

	// std:: cout << "dpsi = " << dpsi << '\n'; 
    // std::cout << "P   " << dpsi << ' ' << Q  << ' ' <<  gs  << ' ' << Assim.a << '\n';
    // std::cout << "P1  " << dpsi << ' ' << Q1 << ' ' <<  gs1 << ' ' << Assim1.a << '\n';

	return (P1-P)/1e-6;

//	double A,B;
//	if (Assim.isVcmaxLimited){
//		A = vcmax; B = Km;
//	}
//	else{
//		double phi0iabs = par_photosynth.phi0 * par_photosynth.Iabs;
//		double jj = 4*phi0iabs/jmax;
//		double jlim = phi0iabs / sqrt(1+ jj*jj);

//		A = jlim; B = 2*gstar;
//	}

//	double x = Assim.ci;

//	double dgs_dchi = (1/ca) * (A*B* (ca - gstar) - br*vcmax*(B + ca*x)*(B + ca*x) + A*ca*(gstar - 2*gstar*x + ca*x*x)) / ((1-x)*(1-x)*(B + ca*x)*(B + ca*x));
//	double dgs_ddpsi = calc_gsprime(dpsi, psi_soil, par_plant, par_env);
//	
//	double dP_ddpsi = dgs_ddpsi * ca * ((1-x) - gs/dgs_dchi) - 2*y*dpsi;

//	return dP_ddpsi;
}

inline double calc_dpsi_bound_inst(double psi_soil, ParPlant par_plant, ParEnv par_env, ParPhotosynth par_photosynth, ParCost par_cost){
  double bound = 100;
  // If using PM, find max dpsi from max possible transpiration 
  if (par_env.et_method == ET_PM){
    double ga = calc_g_aero(par_plant.h_canopy, par_env.v_wind, par_plant.h_wind_measurement);
    double Qmax = calc_max_transpiration_pm(ga, par_env);
    double max_dpsi = calc_dpsi_from_sapflux(Qmax, psi_soil, par_plant, par_env);
    bound = fmin(max_dpsi, bound);
  }

  return bound;
}


} // phydro

#endif
