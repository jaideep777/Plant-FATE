#ifndef PHYDRO_PM_ET_H
#define PHYDRO_PM_ET_H

#include <stdexcept>
#include <cmath>

#include "environment.h"

namespace phydro{



// Aerodynamic conductance [m s-1]
// To convert to mol m-2 s-1, see this: https://rdrr.io/cran/bigleaf/man/ms.to.mol.html (but not convincing)
// Refs: 
//    Eq 13 in Leuning et al (2008). https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2007WR006562
//    Eq 7 in Zhang et al (2008): https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017JD027025
//    Box 4 in https://www.fao.org/3/x0490e/x0490e06.htm 
// h_canopy    canopy height [m]
// v_wind      wind speed [m s-1]
// z_measurement    height of wind speed measurement [m]
inline double calc_g_aero(double h_canopy, double v_wind, double z_measurement){
	double k_karman = 0.41;        // von Karman's constant [-]
	double d = h_canopy*2.0/3.0;   // zero-plane displacement height [m]
	double z_om = 0.123*h_canopy;  // roughness lenghts governing transfer of water and momentum [m]
	double z_ov = 0.1*z_om;
	
	double g_aero = (k_karman*k_karman * v_wind) / ( log((z_measurement-d)/z_om)*log((z_measurement-d)/z_ov) );
	return g_aero;
}


// multiplier to convert:
//   stomatal conductance to CO2 [mol m-2 s-1] ----> stomatal conductance to water [m s-1]
inline double gs_conv(double tc, double patm){
	double R = 8.31446261815324; // Universal gas constant [J mol-1 K-1]

	// PV = nRT, so RT/P = V/n ~ [m3 mol-1]
	// converter from [mol m-2 s-1] to [m s-1]
	return 1.6 * R * (tc+273.16) / patm;
}


// Calculate PML transpiration [mol m-2 s-1]
// gs   Stomatal conductance to CO2 [mol m-2 s-1]
// ga   Aerodynamic conductance [m s-1]
// Rn   Absorbed net radiation [W m-2]
inline double calc_transpiration_pm(double gs, double ga, ParEnv par_env){

	double gw = gs * gs_conv(par_env.tc, par_env.patm);  // gw in [m s-1]

	double latent_energy = (par_env.epsilon*par_env.Rn + (par_env.rho*par_env.cp/par_env.gamma)*ga*par_env.vpd) / (par_env.epsilon + 1 + ga/gw); // latent energy W m-2 
	double trans = latent_energy * (55.5 / par_env.lv); // W m-2 ---> mol m-2 s-1
	return trans;
}


// Calculate maximum possible PML transpiration for a given ga, calculated by setting gs = inf, [mol m-2 s-1]
// ga   Aerodynamic conductance [m s-1]
// Rn   Absorbed net radiation [W m-2]
inline double calc_max_transpiration_pm(double ga, ParEnv par_env){

	double latent_energy = (par_env.epsilon*par_env.Rn + (par_env.rho*par_env.cp/par_env.gamma)*ga*par_env.vpd) / (par_env.epsilon + 1); // latent energy W m-2 
	double trans = latent_energy * (55.5 / par_env.lv); // W m-2 ---> mol m-2 s-1
	return trans;
}


// Calculate PML stomatal conductance to CO2 [mol m-2 s-1]
// Q    Sap flux [mol m-2 s-1]
// ga   Aerodynamic conductance [m s-1]
// Rn   Absorbed net radiation [W m-2]
inline double calc_gs_pm(double Q, double ga, ParEnv par_env){

	double Q_energy = Q * (par_env.lv / 55.5);

	double den = par_env.epsilon*par_env.Rn + (par_env.rho*par_env.cp/par_env.gamma)*ga*par_env.vpd - (1+par_env.epsilon)*Q_energy; 
	//den = fmax(den, 0);

	double gw = ga * Q_energy / den; // stomatal conductance to water [m s-1]

	double gs = gw / gs_conv(par_env.tc, par_env.patm); // stomatal conductance to CO2 [mol m-2 s-1]

	return gs;
}


// Calculate derivative of transpiration wrt stomatal conductance to CO2 [unitless] - analytical version
inline double calc_dE_dgs_pm(double gs, double ga, ParEnv par_env){

	double gw = gs * gs_conv(par_env.tc, par_env.patm);  // [m s-1]

	double num = ga * (par_env.epsilon*par_env.Rn + (par_env.rho*par_env.cp/par_env.gamma)*ga*par_env.vpd);
	double den = par_env.epsilon*gw + gw + ga;

	double d_le_dgw = (num/den/den); // derivative of latent energy wrt stomatal conductance for water in m s-1

	return d_le_dgw * (55.5 / par_env.lv) * gs_conv(par_env.tc, par_env.patm);
}


// Calculate derivative of transpiration wrt stomatal conductance to CO2 [unitless] - numerical version
inline double calc_dE_dgs_pm_num(double gs, double ga, ParEnv par_env){
	double E = calc_transpiration_pm(gs, ga, par_env);
	double E_plus = calc_transpiration_pm(gs+1e-6, ga, par_env);

	return (E_plus-E)/1e-6;
}


} // namespace phydro


#endif


