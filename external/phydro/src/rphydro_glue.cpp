// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       R Interface
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef USINGRCPP

#include <Rcpp.h>
using namespace Rcpp;

#include "phydro.h"

namespace phydro{

// -------------------------------------------------------------
// Utilities to convert strings to enums
// -------------------------------------------------------------

// Conversion function for GsMethod
GsMethod stringToGsMethod(const std::string& methodStr) {
    static const std::unordered_map<std::string, GsMethod> GsMethod_map {
        {"GS_IGF", GS_IGF},
        {"GS_QNG", GS_QNG},
        {"GS_APX", GS_APX},
        {"GS_APX2", GS_APX2}
    };

    return GsMethod_map.at(methodStr);
}

// Conversion function for ETMethod
ETMethod stringToETMethod(const std::string& methodStr) {
    static const std::unordered_map<std::string, ETMethod> ETMethod_map {
        {"ET_DIFFUSION", ET_DIFFUSION},
        {"ET_PM", ET_PM}
    };

    return ETMethod_map.at(methodStr);
}

// Conversion function for FtempVcmaxJmaxMethod
FtempVcmaxJmaxMethod stringToFtempVcmaxJmaxMethod(const std::string& methodStr) {
    static const std::unordered_map<std::string, FtempVcmaxJmaxMethod> FtempVcmaxJmaxMethod_map {
        {"FV_kattge07", FV_kattge07},
        {"FV_kumarathunge19", FV_kumarathunge19},
        {"FV_leuning02", FV_leuning02}
    };

    return FtempVcmaxJmaxMethod_map.at(methodStr);
}

// Conversion function for FtempRdMethod
FtempRdMethod stringToFtempRdMethod(const std::string& methodStr) {
    static const std::unordered_map<std::string, FtempRdMethod> FtempRdMethod_map {
        {"FR_heskel16", FR_heskel16},
        {"FR_arrhenius", FR_arrhenius},
        {"FR_q10", FR_q10}
    };

    return FtempRdMethod_map.at(methodStr);
}

// Conversion function for FtempBrMethod
FtempBrMethod stringToFtempBrMethod(const std::string& methodStr) {
    static const std::unordered_map<std::string, FtempBrMethod> FtempBrMethod_map {
        {"FB_atkin15", FB_atkin15},
        {"FB_kumarathunge19", FB_kumarathunge19}
    };

    return FtempBrMethod_map.at(methodStr);
}

// -------------------------------------------------------------
//    Utility to convert Phydro result to R list
// -------------------------------------------------------------
inline Rcpp::List PHydroResult_to_List(const phydro::PHydroResult& res){
	return Rcpp::List::create(
	           Named("a") = res.a,
	           Named("e") = res.e,
	           Named("gs") = res.gs,
	           Named("ci") = res.ci,
	           Named("chi") = res.chi,
	           Named("vcmax") = res.vcmax,
	           Named("jmax") = res.jmax,
	           Named("dpsi") = res.dpsi,
//	           Named("psi_l") = res.psi_l,
	           Named("mc") = res.mc,
	           Named("mj") = res.mj,
	           Named("gammastar") = res.gammastar,
	           Named("kmm") = res.kmm,
	           Named("vcmax25") = res.vcmax25,
	           Named("jmax25") = res.jmax25,
	           Named("rd") = res.rd,
	           Named("isVcmaxLimited") = res.isVcmaxLimited,
	           Named("ac") = res.ac,
	           Named("aj") = res.aj,
			   Named("le") = res.le,
			   Named("le_s_wet") = res.le_s_wet
	       );
}

// -------------------------------------------------------------
//   Utility to convert control options list to struct
// -------------------------------------------------------------
ParControl listToParControl(Rcpp::List options){
	ParControl par_control;
	par_control.gs_method       = stringToGsMethod(options["gs_method"]);
	par_control.et_method       = stringToETMethod(options["et_method"]);
	par_control.ftemp_vj_method = stringToFtempVcmaxJmaxMethod(options["ftemp_vj_method"]);
	par_control.ftemp_rd_method = stringToFtempRdMethod(options["ftemp_rd_method"]);
	par_control.ftemp_br_method = stringToFtempBrMethod(options["ftemp_br_method"]);
	par_control.scale_alpha     = options["scale_alpha"];
	return par_control;
}


ParPlant listToParPlant(Rcpp::List l){
	ParPlant par_plant_cpp(l["conductivity"], l["psi50"], l["b"]);

	if (l.containsElementNamed("tchome")) par_plant_cpp.tchome = l["tchome"];
	if (l.containsElementNamed("h_canopy")) par_plant_cpp.h_canopy = l["h_canopy"];
	if (l.containsElementNamed("h_wind_measurement")) par_plant_cpp.h_wind_measurement = l["h_wind_measurement"];
	return par_plant_cpp;
}


// Wrappers for temperature response functions
inline double r_calc_ftemp_inst_jmax(double tcleaf, double tcgrowth, double tchome, double tcref, std::string method_ftemp) {
	return calc_ftemp_inst_jmax(tcleaf, tcgrowth, tchome, tcref, stringToFtempVcmaxJmaxMethod(method_ftemp));
}

inline double r_calc_ftemp_inst_vcmax(double tcleaf, double tcgrowth, double tcref, std::string method_ftemp) {
	return calc_ftemp_inst_vcmax(tcleaf, tcgrowth, tcref, stringToFtempVcmaxJmaxMethod(method_ftemp));
}


// -------------------------------------------------------------
//  Wrappers for Phydro calls taking R lists instead of parameter objects
// -------------------------------------------------------------
inline Rcpp::List rphydro_analytical(double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, Rcpp::List par_plant, Rcpp::List par_cost, Rcpp::List options){
	ParControl par_control = listToParControl(options);
	ParPlant par_plant_cpp = listToParPlant(par_plant);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);

	return PHydroResult_to_List(phydro_analytical(tc, tg, ppfd, netrad, vpd, co2, pa, fapar, kphio, psi_soil, rdark, vwind, par_plant_cpp, par_cost_cpp, par_control));
}

inline Rcpp::List rphydro_instantaneous_analytical(double vcmax25, double jmax25, double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, Rcpp::List par_plant, Rcpp::List par_cost, Rcpp::List options){
	ParControl par_control = listToParControl(options);
	ParPlant par_plant_cpp = listToParPlant(par_plant);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);

	return PHydroResult_to_List(phydro_instantaneous_analytical(vcmax25, jmax25, tc, tg, ppfd, netrad, vpd, co2, pa, fapar, kphio, psi_soil, rdark, vwind, par_plant_cpp, par_cost_cpp, par_control));
}


#ifndef PHYDRO_ANALYTICAL_ONLY

inline Rcpp::List rphydro_numerical(double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, Rcpp::List par_plant, Rcpp::List par_cost, Rcpp::List options){
	ParControl par_control = listToParControl(options);
	ParPlant par_plant_cpp = listToParPlant(par_plant);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);

	return PHydroResult_to_List(phydro_numerical(tc, tg, ppfd, netrad, vpd, co2, pa, fapar, kphio, psi_soil, rdark, vwind, par_plant_cpp, par_cost_cpp, par_control));
}

inline Rcpp::List rphydro_instantaneous_numerical(double vcmax25, double jmax25, double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, Rcpp::List par_plant, Rcpp::List par_cost, Rcpp::List options){
	ParControl par_control = listToParControl(options);
	ParPlant par_plant_cpp = listToParPlant(par_plant);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);

	return PHydroResult_to_List(phydro_instantaneous_numerical(vcmax25, jmax25, tc, tg, ppfd, netrad, vpd, co2, pa, fapar, kphio, psi_soil, rdark, vwind, par_plant_cpp, par_cost_cpp, par_control));
}

#endif


} // phydro namespace


RCPP_EXPOSED_CLASS_NODECL(phydro::ParEnv);
RCPP_EXPOSED_CLASS_NODECL(phydro::ParPlant);

using namespace phydro;

// -------------------------------------------------------------
//   R Interface
// -------------------------------------------------------------
RCPP_MODULE(phydro_module) {
	class_ <phydro::ParEnv>("ParEnv")
		.constructor<double,double,double,double,double>()
		.constructor<double,double,double,double>()
		.method("print", &ParEnv::print)
		.method("calc_temp_dependencies", &ParEnv::calc_temp_dependencies)

		.field("tc", &ParEnv::tc)
		.field("vpd", &ParEnv::vpd)
		.field("patm", &ParEnv::patm)
		.field("Rn", &ParEnv::Rn)
		.field("v_wind", &ParEnv::v_wind)
	;

	class_ <phydro::ParPlant>("ParPlant")
		.constructor<double,double,double,double,double,double>()
		.constructor<double,double,double>()
		.method("print", &ParPlant::print)

		.field("conductivity", &ParPlant::conductivity)
		.field("psi50", &ParPlant::psi50)
		.field("b", &ParPlant::b)
		.field("tchome", &ParPlant::tchome)
		.field("h_canopy", &ParPlant::h_canopy)
		.field("h_wind_measurement", &ParPlant::h_wind_measurement)
	;

	function("listToParPlant", &listToParPlant);

	// Temperature dependencies
	function("calc_kmm", &calc_kmm);
	function("calc_patm", &calc_patm);
	function("calc_gammastar", &calc_gammastar);
	function("calc_density_h2o", &calc_density_h2o);
	function("calc_viscosity_h2o", &calc_viscosity_h2o);
	function("calc_ftemp_kphio", &calc_ftemp_kphio);
	function("r_calc_ftemp_inst_vcmax", &r_calc_ftemp_inst_vcmax);
	function("r_calc_ftemp_inst_jmax", &r_calc_ftemp_inst_jmax);
	// function("calc_ftemp_inst_rd", &calc_ftemp_inst_rd);

	// pm functions
	function("calc_esat", &calc_esat);
	function("calc_g_aero", &calc_g_aero);
	function("calc_density_air", &calc_density_air);
	function("calc_enthalpy_vap", &calc_enthalpy_vap);
	function("calc_cp_moist_air", &calc_cp_moist_air);
	function("calc_psychro", &calc_psychro);
	function("calc_sat_slope", &calc_sat_slope);
	function("calc_transpiration_pm", &calc_transpiration_pm);
	function("calc_max_transpiration_pm", &calc_max_transpiration_pm);
	function("calc_gs_pm", &calc_gs_pm);
	function("calc_dE_dgs_pm", &calc_dE_dgs_pm);
	function("calc_dE_dgs_pm_num", &calc_dE_dgs_pm_num);


	// Phydro core
	function("rphydro_analytical", &rphydro_analytical);
	function("rphydro_instantaneous_analytical", &rphydro_instantaneous_analytical);

#ifndef PHYDRO_ANALYTICAL_ONLY
	function("rphydro_numerical", &rphydro_numerical);
	function("rphydro_instantaneous_numerical", &rphydro_instantaneous_numerical);
#endif

}

#endif

