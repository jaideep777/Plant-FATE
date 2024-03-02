#ifndef PHYDRO_ANALYTICAL_SOLVER_H
#define PHYDRO_ANALYTICAL_SOLVER_H

#include "pn_zero.h"

#include "hyd_transpiration.h"
#include "hyd_photosynthesis.h"


namespace phydro{

inline double calc_J(double gs, double x, ParPhotosynth par_photosynth){
	double g = par_photosynth.gammastar / par_photosynth.ca;
	double k = par_photosynth.kmm / par_photosynth.ca;
	double ca = par_photosynth.ca / par_photosynth.patm*1e6;
	double d = par_photosynth.delta;
	return 4*gs*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k));
}


inline double calc_jmax_from_J(double J, ParPhotosynth par_photosynth){
	double p = 4*par_photosynth.phi0 * par_photosynth.Iabs;
	double pj = p/J;
	return p/sqrt(pj*pj-1);
}


inline double calc_djmax_dJ(double J, ParPhotosynth par_photosynth){
  double p = 4*par_photosynth.phi0 * par_photosynth.Iabs;
  double sq = sqrt(p*p-J*J);
  double psq = p/sq;
  return psq*psq*psq;
  //return 4*p^3/ajmax^3/((p/ajmax)^2-1)^(3/2);
}


inline double calc_dJ_dchi(double gs, double x, ParPhotosynth par_photosynth){
  double g = par_photosynth.gammastar / par_photosynth.ca;
  double k = par_photosynth.kmm / par_photosynth.ca;
  double ca = par_photosynth.ca / par_photosynth.patm*1e6;
  double d = par_photosynth.delta;
  // gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) + 2*g^2 + g*(2*x - 3) - x^2)/(d*(k + x) + g - x)^2)
  double D = d*(k + x) + g - x;
  return 4*gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x*x) - ((x-g)*(x-g)+3*g*(1-g)))/(D*D));
  // gs*ca*(3*(g-1)*g/(g-x)^2 - 1)
}


inline double calc_dJ_ddpsi(double gsprime, double x, ParPhotosynth par_photosynth){
  double g = par_photosynth.gammastar/par_photosynth.ca;
  double k = par_photosynth.kmm/par_photosynth.ca;
  double ca = par_photosynth.ca/par_photosynth.patm*1e6;
  double d = par_photosynth.delta;
  return 4*gsprime*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k));
}


inline double calc_x_from_dpsi(double dpsi, double gsprime, ParPhotosynth par_photosynth, ParCost par_cost){
  double gstar = par_photosynth.gammastar/par_photosynth.patm*1e6;
  double Km = par_photosynth.kmm/par_photosynth.patm*1e6;
  double ca = par_photosynth.ca/par_photosynth.patm*1e6;
  double br = par_photosynth.delta;
  double y = par_cost.gamma;
  
  double ca2 = ca*ca;
  double x = (-2*ca*dpsi*(gstar + br*Km)*y + 
     ca2*((3 - 2*br)*gstar + br*Km)*gsprime + 
     -sqrt(2)*sqrt(
       ca2*dpsi*((-3 + 2*br)*gstar - br*Km)*((-1 + br)*ca + gstar + 
                                                br*Km)*y*
         (-2*dpsi*y + (ca + 2*gstar)*
            gsprime)))/
    (ca2*(2*(-1 + br)*dpsi*y + ((3 - 2*br)*gstar + br*Km)*
             gsprime));
  
  if (x < (gstar + br*Km)/(ca - br*ca)) x = (gstar + br*Km)/(ca - br*ca)+1e-12;
  return x;
}


// inline double calc_delta_from_dpsi(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env, ParPhotosynth par_photosynth, ParCost par_cost){
//   double gstar = par_photosynth.gammastar/par_photosynth.patm*1e6;
//   double Km = par_photosynth.kmm/par_photosynth.patm*1e6;
//   double ca = par_photosynth.ca/par_photosynth.patm*1e6;
//   double br = par_photosynth.delta;
//   double y = par_cost.gamma;
   
//   double gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env);
  
//   double delt = (-2*dpsi*y + (ca + 2*gstar)*gsprime);
//   return delt;
// }


struct DPsiBounds{
	double exact;
	double approx_O2;
	double Iabs_bound;
	DPsiBounds(double e, double a, double i) : exact(e), approx_O2(a), Iabs_bound(i) {}
};

inline DPsiBounds calc_dpsi_bound(double psi_soil, ParPlant par_plant, ParEnv par_env, ParPhotosynth par_photosynth, ParCost par_cost){
  double gstar = par_photosynth.gammastar/par_photosynth.patm*1e6;
  double ca = par_photosynth.ca/par_photosynth.patm*1e6;
  double y = par_cost.gamma;
  double K = scale_conductivity(par_plant.conductivity, par_env)/(1.6*par_env.vpd/par_env.patm);
  //double K = K/(1.6*par_env.vpd/par_env.patm);
  double Pox = P(psi_soil, par_plant.psi50, par_plant.b);
  double Ppox = Pprime(psi_soil, par_plant.psi50, par_plant.b);
  double Pppox = Pprimeprime(psi_soil, par_plant.psi50, par_plant.b);
  
  auto f1 = [&](double dpsi){
    double Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env);
    double gs = calc_gs_from_Q(Q, psi_soil, par_plant, par_env);
    double gsprime = calc_gsprime(dpsi, gs, psi_soil, par_plant, par_env);
    double x = calc_x_from_dpsi(dpsi,gsprime, par_photosynth, par_cost);
    double J=calc_J(gs, x, par_photosynth)-4*par_photosynth.phi0*par_photosynth.Iabs;
    return J;
  };
  
  double a = (ca + 2*gstar)*K*Pppox*4/8;
  double b = -(2*y + (ca + 2*gstar)*K*Ppox);
  double c = (ca + 2*gstar)*K*Pox;
  double del = b*b-4*a*c;

  auto f2 = [&](double dpsi){
    double gsprime = calc_gsprime_from_dpsi(dpsi, psi_soil, par_plant, par_env);
    return (-2*dpsi*y + (ca + 2*gstar)*gsprime);
  };

  double approx_O2 = (-b-sqrt(del))/(2*a);
  double exact = pn::zero(0,10, f2, 1e-6).root;

  double use_bound = exact;
	
  //# cat(psi_soil, ":", exact, " ", approx_O2, " ", use_bound, "\n");
  double Iabs_bound = pn::zero(use_bound*0.001, use_bound*0.99, f1, 1e-6).root;
  
  // If using PM, find max dpsi from max possible transpiration 
  if (par_env.et_method == ET_PM){
    double ga = calc_g_aero(par_plant.h_canopy, par_env.v_wind, par_plant.h_wind_measurement);
    double Qmax = calc_max_transpiration_pm(ga, par_env);
    double max_dpsi = calc_dpsi_from_sapflux(Qmax, psi_soil, par_plant, par_env);
    Iabs_bound = fmin(max_dpsi, Iabs_bound);
  }

  //# dpsi=seq(exact*0.001,exact*0.99, length.out=200);
  //# plot(y=sapply(X = dpsi, FUN = f1), x=dpsi, type="l");
  
  return DPsiBounds(exact, approx_O2, Iabs_bound);
}


struct DFDX{
	double dPdx;
	double J;
	double djmax_dJ;
	double dJ_dchi;
	DFDX(double _dpdx, double _J, double _djmaxdJ, double _dJdx): dPdx(_dpdx), J(_J), djmax_dJ(_djmaxdJ), dJ_dchi(_dJdx) {}
};


inline DFDX dFdx(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env, ParPhotosynth par_photosynth, ParCost par_cost){
  
  double Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env);
  double gs = calc_gs_from_Q(Q, psi_soil, par_plant, par_env);
  double gsprime = calc_gsprime(dpsi, gs, psi_soil, par_plant, par_env);
  
  double X =  calc_x_from_dpsi(dpsi, gsprime, par_photosynth, par_cost);
  
  double J = calc_J(gs, X, par_photosynth);
  
  double ca = par_photosynth.ca / par_photosynth.patm*1e6;
  double g = par_photosynth.gammastar / par_photosynth.ca;
  
  double djmax_dJ = calc_djmax_dJ(J, par_photosynth);
  double dJ_dchi  = calc_dJ_dchi(gs, X, par_photosynth);
  
  double dP_dx = -gs*ca - par_cost.alpha * djmax_dJ * dJ_dchi;
  
  //# dP_ddpsi = gsprime*ca*(1-X) - par_cost$alpha * calc_djmax_dAjmax(ajmax, par_photosynth) * calc_dAjmax_ddpsi(gsprime, X, par_photosynth) - 2*par_cost$gamma*dpsi #/par_plant$psi50^2
  //# # cat(c(dP_dx, dP_ddpsi), "\n")
  //# c(dP_dx, dP_ddpsi)
  return DFDX(dP_dx, J, djmax_dJ, dJ_dchi);
}


inline double chi_jmax_limited(ParPhotosynth par_photosynth, ParCost par_cost){
  double g = par_photosynth.gammastar/par_photosynth.ca;
  double k = par_photosynth.kmm/par_photosynth.ca;
  double b = par_photosynth.delta;
  double a = par_cost.alpha;
  
  //#(g*(1-4*a) + 2*sqrt(3)*sqrt(a*(1-4*a)*(1-g)*g))/(1-4*a)
  //# (2*sqrt(-a*(4*a + b - 1)*(2*b^2*g*k + 2*b^2*g + b^2*(-k^2) - b^2*k + 2*b*g^2 - 4*b*g*k - 5*b*g + b*k - 3*g^2 + 3*g)) - 4*a*b*k - 4*a*g + b^2*(-k) - b*g + b*k + g)/((b - 1)*(4*a + b - 1))
  return (2*sqrt(-a*(4*a + b - 1)*(-3*g + 2*b*g - b*k)*(-1 + b + g + b*k)) - (4*a + b - 1)*(b*k + g))/((b - 1)*(4*a + b - 1));
}



} // phydro

#endif
