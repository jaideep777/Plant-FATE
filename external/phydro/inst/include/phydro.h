#ifndef PHYDRO_PHYDRO_H
#define PHYDRO_PHYDRO_H


#include "hyd_analytical_solver.h"
#include "hyd_instantaneous_solver_analytical.h"

namespace phydro{

struct PHydroResult{
	double a;
	double e;
	double gs;
	double ci;
	double chi;
	double vcmax;
	double jmax;
	double dpsi;
	double psi_l;
	double nfnct;
	double niter;
	double mc;
	double mj;
	double gammastar;
	double kmm;
	double vcmax25;
	double jmax25;
	double rd;
	bool   isVcmaxLimited;
	double ac;
	double aj;
	double le;
	double le_s_wet;
};


class ParControl{
	public:
	GsMethod             gs_method       = GS_IGF;
	ETMethod             et_method       = ET_DIFFUSION;
	FtempVcmaxJmaxMethod ftemp_vj_method = FV_kumarathunge19;
	FtempRdMethod        ftemp_rd_method = FR_heskel16; 
	FtempBrMethod        ftemp_br_method = FB_atkin15;
	bool                 scale_alpha     = false;
};


inline PHydroResult phydro_analytical(double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, ParPlant par_plant, ParCost par_cost = ParCost(0.1,1), ParControl par_control = ParControl()){
	
	// double pa = calc_patm(elv);

	ParPhotosynth par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark, tg, par_plant.tchome, par_control.ftemp_vj_method, par_control.ftemp_rd_method, par_control.ftemp_br_method);
	ParEnv        par_env(tc, pa, vpd, netrad, vwind);
	if (par_control.scale_alpha) par_cost.alpha /= par_photosynth.fT_jmax; // convert alpha from cost of jmax to cost of jmax25
	par_env.gs_method = par_control.gs_method;
	par_env.et_method = par_control.et_method;

	auto     bounds = calc_dpsi_bound(psi_soil, par_plant, par_env, par_photosynth, par_cost);
	auto   dpsi_opt = pn::zero(bounds.Iabs_bound*0.001, bounds.Iabs_bound*0.999, [&](double dpsi){return dFdx(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost).dPdx;}, 1e-6);
	double        e = calc_sapflux(dpsi_opt.root, psi_soil, par_plant, par_env);
	double       gs = calc_gs_from_Q(e, psi_soil, par_plant, par_env);
    double  gsprime = calc_gsprime(dpsi_opt.root, gs, psi_soil, par_plant, par_env);
    double        x = calc_x_from_dpsi(dpsi_opt.root,gsprime,par_photosynth, par_cost);
	double        J = calc_J(gs, x, par_photosynth);	
	double     jmax = calc_jmax_from_J(J, par_photosynth); 
	double    vcmax = (J/4.0)*(x*par_photosynth.ca + par_photosynth.kmm)/(x*par_photosynth.ca + 2*par_photosynth.gammastar);
	double        a = gs*(par_photosynth.ca/par_photosynth.patm*1e6)*(1-x);

	PHydroResult res;
	res.a = a;
	res.e = e; //1.6*gs*vpd/par_env.patm;
	res.ci = x*par_photosynth.ca;
	res.gs = gs;
	res.chi = x;
	res.vcmax = vcmax;
	res.jmax = jmax;
	res.dpsi = dpsi_opt.root;
	res.psi_l = psi_soil - dpsi_opt.root;
	res.nfnct = dpsi_opt.nfnct;
	res.mc = (x*par_photosynth.ca - par_photosynth.gammastar) / (x*par_photosynth.ca + par_photosynth.kmm);
	res.mj = (x*par_photosynth.ca - par_photosynth.gammastar) / (x*par_photosynth.ca + 2*par_photosynth.gammastar);
	res.gammastar = par_photosynth.gammastar;
	res.kmm = par_photosynth.kmm;
	res.vcmax25 = vcmax / par_photosynth.fT_vcmax; // calc_ftemp_vcmax_bernacchi(tc);
	res.jmax25 = jmax / par_photosynth.fT_jmax;
	res.rd = vcmax * par_photosynth.delta;
	res.isVcmaxLimited = 0.5; // coordinated
	res.ac = a;
	res.aj = a;
	res.le = e * 0.018015 * par_env.lv; // e [mol m-2 s-1] * 0.01815 [kg mol-1] * lv [J kg-1] = J m-2 s-1 = W m-2
	res.le_s_wet = (1-fapar)*netrad*(par_env.epsilon/(1+par_env.epsilon));

	return res;

}


inline PHydroResult phydro_instantaneous_analytical(double vcmax25, double jmax25, double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, ParPlant par_plant, ParCost par_cost = ParCost(0.1,1), ParControl par_control = ParControl()){
	
	// double pa = calc_patm(elv);

	ParPhotosynth par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark, tg, par_plant.tchome, par_control.ftemp_vj_method, par_control.ftemp_rd_method, par_control.ftemp_br_method);
	ParEnv        par_env(tc, pa, vpd, netrad, vwind);
	if (par_control.scale_alpha) par_cost.alpha /= par_photosynth.fT_jmax; // convert alpha from cost of jmax to cost of jmax25
	par_env.gs_method = par_control.gs_method;
	par_env.et_method = par_control.et_method;

	// par_photosynth.print();

	double vcmax = vcmax25 * par_photosynth.fT_vcmax;
	double jmax  =  jmax25 * par_photosynth.fT_jmax;

	double  bound = calc_dpsi_bound_inst(psi_soil, par_plant, par_env, par_photosynth, par_cost);
	auto dpsi_opt = pn::zero(0, 0.99*bound, [&](double dpsi){return calc_dP_ddpsi(dpsi, vcmax, jmax, psi_soil, par_plant, par_env, par_photosynth, par_cost);}, 1e-6);
	double      e = calc_sapflux(dpsi_opt.root, psi_soil, par_plant, par_env);
	double     gs = calc_gs_from_Q(e, psi_soil, par_plant, par_env);
	auto        A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth); 	

	PHydroResult res;
	res.a = A.a;
	res.e = e; //1.6*gs*vpd/par_env.patm;
	res.ci = A.ci;
	res.gs = gs;
	res.chi = A.ci/par_photosynth.ca;
	res.vcmax = vcmax;
	res.jmax = jmax;
	res.dpsi = dpsi_opt.root;
	res.psi_l = psi_soil - dpsi_opt.root;
	res.mc = (A.ci - par_photosynth.gammastar) / (A.ci + par_photosynth.kmm);
	res.mj = (A.ci - par_photosynth.gammastar) / (A.ci + 2*par_photosynth.gammastar);
	res.gammastar = par_photosynth.gammastar;
	res.kmm = par_photosynth.kmm;
	res.vcmax25 = vcmax25;
	res.jmax25 = jmax25;
	res.rd = vcmax * par_photosynth.delta;
	res.isVcmaxLimited = A.isVcmaxLimited;
	res.ac = calc_assim_rubisco_limited(gs, vcmax, par_photosynth).a;
	res.aj = calc_assim_light_limited(gs, jmax, par_photosynth).a; 	
	res.le = e * 0.018015 * par_env.lv; // e [mol m-2 s-1] * 0.01815 [kg mol-1] * lv [J kg-1] = J m-2 s-1 = W m-2
	res.le_s_wet = (1-fapar)*netrad*(par_env.epsilon/(1+par_env.epsilon));

	return res;

}

} // phydro


#ifndef PHYDRO_ANALYTICAL_ONLY

#include "hyd_numerical_solver.h"
#include "hyd_instantaneous_solver_numerical.h"


namespace phydro{

inline PHydroResult phydro_numerical(double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, ParPlant par_plant, ParCost par_cost = ParCost(0.1,1), ParControl par_control = ParControl()){
	
	// double pa = calc_patm(elv);

	ParPhotosynth par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark, tg, par_plant.tchome, par_control.ftemp_vj_method, par_control.ftemp_rd_method, par_control.ftemp_br_method);
	ParEnv        par_env(tc, pa, vpd, netrad, vwind);
	if (par_control.scale_alpha) par_cost.alpha /= par_photosynth.fT_jmax; // convert alpha from cost of jmax to cost of jmax25
	par_env.gs_method = par_control.gs_method;
	par_env.et_method = par_control.et_method;

	auto     opt = optimize_midterm_multi(psi_soil, par_cost, par_photosynth, par_plant, par_env);
	double     e = calc_sapflux(opt.dpsi, psi_soil, par_plant, par_env);
	double    gs = calc_gs_from_Q(e, psi_soil, par_plant, par_env);
	auto      aj = calc_assim_light_limited(gs, opt.jmax, par_photosynth); 	
	double vcmax = vcmax_coordinated_numerical(aj.a, aj.ci, par_photosynth);

	PHydroResult res;
	res.a = aj.a;
	res.e = e; //1.6*gs*vpd/par_env.patm;
	res.ci = aj.ci;
	res.gs = gs;
	res.chi = aj.ci/par_photosynth.ca;
	res.vcmax = vcmax;
	res.jmax = opt.jmax;
	res.dpsi = opt.dpsi;
	res.psi_l = psi_soil - opt.dpsi;
	res.nfnct = 0;
	res.mc = (aj.ci - par_photosynth.gammastar) / (aj.ci + par_photosynth.kmm);
	res.mj = (aj.ci - par_photosynth.gammastar) / (aj.ci + 2*par_photosynth.gammastar);
	res.gammastar = par_photosynth.gammastar;
	res.kmm = par_photosynth.kmm;
	res.vcmax25 = vcmax / par_photosynth.fT_vcmax; // calc_ftemp_vcmax_bernacchi(tc);
	res.jmax25 = opt.jmax / par_photosynth.fT_jmax;
	res.rd = vcmax * par_photosynth.delta;
	res.isVcmaxLimited = 0.5; // coordinated
	res.ac = aj.a;
	res.aj = aj.a;
	res.le = e * 0.018015 * par_env.lv; // e [mol m-2 s-1] * 0.01815 [kg mol-1] * lv [J kg-1] = J m-2 s-1 = W m-2
	res.le_s_wet = (1-fapar)*netrad*(par_env.epsilon/(1+par_env.epsilon));

	return res;

}


inline PHydroResult phydro_instantaneous_numerical(double vcmax25, double jmax25, double tc, double tg, double ppfd, double netrad, double vpd, double co2, double pa, double fapar, double kphio, double psi_soil, double rdark, double vwind, ParPlant par_plant, ParCost par_cost = ParCost(0.1,1), ParControl par_control = ParControl()){
	
	// double pa = calc_patm(elv);

	ParPhotosynth par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark, tg, par_plant.tchome, par_control.ftemp_vj_method, par_control.ftemp_rd_method, par_control.ftemp_br_method);
	ParEnv        par_env(tc, pa, vpd, netrad, vwind);
	if (par_control.scale_alpha) par_cost.alpha /= par_photosynth.fT_jmax; // convert alpha from cost of jmax to cost of jmax25
	par_env.gs_method = par_control.gs_method;
	par_env.et_method = par_control.et_method;

	double vcmax = vcmax25 * par_photosynth.fT_vcmax;
	double jmax  =  jmax25 * par_photosynth.fT_jmax;

	double  dpsi = optimize_shortterm_multi(vcmax, jmax, psi_soil, par_cost, par_photosynth, par_plant, par_env);
	double     e = calc_sapflux(dpsi, psi_soil, par_plant, par_env);
	double    gs = calc_gs_from_Q(e, psi_soil, par_plant, par_env);
	auto       A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth); 	

	PHydroResult res;
	res.a = A.a;
	res.e = e; //1.6*gs*vpd/par_env.patm;
	res.ci = A.ci;
	res.gs = gs;
	res.chi = A.ci/par_photosynth.ca;
	res.vcmax = vcmax;
	res.jmax = jmax;
	res.dpsi = dpsi;
	res.psi_l = psi_soil - dpsi;
	res.mc = (A.ci - par_photosynth.gammastar) / (A.ci + par_photosynth.kmm);
	res.mj = (A.ci - par_photosynth.gammastar) / (A.ci + 2*par_photosynth.gammastar);
	res.gammastar = par_photosynth.gammastar;
	res.kmm = par_photosynth.kmm;
	res.vcmax25 = vcmax25;
	res.jmax25 = jmax25;
	res.rd = vcmax * par_photosynth.delta;
	res.isVcmaxLimited = A.isVcmaxLimited;
	res.ac = calc_assim_rubisco_limited(gs, vcmax, par_photosynth).a;
	res.aj = calc_assim_light_limited(gs, jmax, par_photosynth).a; 	
	res.le = e * 0.018015 * par_env.lv; // e [mol m-2 s-1] * 0.01815 [kg mol-1] * lv [J kg-1] = J m-2 s-1 = W m-2
	res.le_s_wet = (1-fapar)*netrad*(par_env.epsilon/(1+par_env.epsilon));

	return res;

}

} // phydro


#endif   // PHYDRO_ANALYTICAL_ONLY


#endif  // include guard
