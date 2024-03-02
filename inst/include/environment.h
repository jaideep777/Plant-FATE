#ifndef PHYDRO_ENVIRONMENT_H
#define PHYDRO_ENVIRONMENT_H

#include <iostream>
#include "temperature_dependencies_physical.h"

namespace phydro{

enum GsMethod{GS_IGF, GS_QNG, GS_APX, GS_APX2};
enum ETMethod{ET_DIFFUSION, ET_PM};

class ParEnv{
	public:
	double tc;     // Temperature [degC]
	double patm;   // Atmospheric pressure [Pa]
	double vpd;    // VPD [Pa]
	double Rn;     // net radiation [W m-2]
	double v_wind; // Wind speed [m s-1]

	double viscosity_water;  // [Pa s]
	double density_water;    // [kg m-3]

	double rho;       // density of air [kg m-3]
	double cp;	      // specific heat capacity of moist air [J kg-1 K-1]
	double gamma;     // psychrometric constant [Pa K-1]
	double epsilon;   // slope of saturation-pressure - temp curve [Pa K-1]
	double lv;        // latent heat of vaporization of water [J kg-1]


	GsMethod gs_method = GS_IGF;
	ETMethod et_method = ET_DIFFUSION;

	inline ParEnv(double _tc, double _patm, double _vpd, double _Rn, double _v_wind){
		tc = _tc;
		vpd = _vpd;
		patm = _patm;
		Rn = _Rn;
		v_wind = _v_wind; 

		calc_temp_dependencies();
	}

	// Note: Using separate constructor instead of default value for _v_wind because Rcpp cannot handle default values
	inline ParEnv(double _tc, double _patm, double _vpd, double _Rn) : ParEnv(_tc, _patm, _vpd, _Rn, 3) { // global average value of v_wind
	}


	inline void calc_temp_dependencies(){
		viscosity_water = calc_viscosity_h2o(tc, patm);
		density_water = calc_density_h2o(tc, patm);

		rho = calc_density_air(tc, patm, vpd, true);
		cp = calc_cp_moist_air(tc);
		gamma = calc_psychro(tc, patm);
		epsilon = calc_sat_slope(tc) / gamma;

		lv = calc_enthalpy_vap(tc);
	}

	inline void print(){
		std::cout << "Env:\n";
		std::cout << "   tc = " << tc << " [degC]\n";     // Temperature [degC]
		std::cout << "   patm = " << patm << " [Pa]\n";   // Atmospheric pressure [Pa]
		std::cout << "   vpd = " << vpd << " [Pa]\n";    // VPD [Pa]
		std::cout << "   Rn = " << Rn << " [w m-2]\n";     // net radiation [w m-2]
		std::cout << "   v_wind = " << v_wind << " [m s-1]\n"; // Wind speed [m s-1]
		std::cout << "   viscosity_water = " << viscosity_water << " [Pa s]\n";  // [Pa s]
		std::cout << "   density_water = " << density_water << " [kg m-3]\n";    // [kg m-3]
		std::cout << "   rho = " << rho << " [kg m-3]\n";       // density of air [kg m-3]
		std::cout << "   cp = " << cp << " [J kg-1 K-1]\n";	      // specific heat capacity of moist air 
		std::cout << "   gamma = " << gamma << " [Pa K-1]\n";     // psychrometric constant
		std::cout << "   epsilon = " << epsilon << " [Pa K-1]\n";   // slope of saturation-pressure - temp curve
		std::cout << "   lv = " << lv << " [J kg-1]\n";        // latent heat of vaporization of water
	}
};


} // namespace phydro

#endif


