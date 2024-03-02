#ifndef PHYDRO_INST_NUMERICAL_SOLVER_H
#define PHYDRO_INST_NUMERICAL_SOLVER_H

#ifdef USINGRCPP
#include <RcppEigen.h>
#else
#include <Eigen/Core>
#endif
 
#include <LBFGSB.h>
#include <LBFGS.h>

#include <iostream>

#include "hyd_transpiration.h"
#include "hyd_photosynthesis.h"

using Eigen::VectorXd;

namespace phydro{

class PHydro_Profit_Inst{
	private:

	int n = 1;

	double psi_soil, vcmax, jmax;

	ParCost       par_cost;
	ParEnv        par_env;
	ParPhotosynth par_photosynth;
	ParPlant      par_plant;

	public:

	inline PHydro_Profit_Inst(double _vcmax, double _jmax, double _psi_soil, ParCost _par_cost, ParPhotosynth _par_photosynth, ParPlant _par_plant, ParEnv _par_env) : 
		psi_soil       ( _psi_soil),
		vcmax          ( _vcmax),
		jmax           ( _jmax),
		par_cost       ( _par_cost),
		par_env        ( _par_env),
		par_photosynth ( _par_photosynth),
		par_plant      ( _par_plant) {
	}

	inline double value(const VectorXd &x) {
		double dpsi = x[0];
		
		double Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env);
		double gs = calc_gs_from_Q(Q, psi_soil, par_plant, par_env);
		auto   A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth);  // min(Ac, Aj) in umol/m2/s
		
		double costs =  par_cost.gamma * dpsi*dpsi;

		double profit = A.a - costs;
		
		//std::cout << "dpsi = " << dpsi << ", profit = " << profit << std::endl;
		return -profit;
	}

    inline double operator()(const VectorXd& x, VectorXd& grad){
		double f = value(x);
		
		for (int i=0; i<n; ++i){
			VectorXd dx = VectorXd::Zero(n);
			double delta = 2.2204e-6;
			dx[i] = delta;
			double fplus = value(x+dx);
			double fminus = value(x-dx);
			grad[i] = (fplus-fminus)/delta/2;
		}

		return f;
	}
};


inline double optimize_shortterm_multi(double vcmax, double jmax, double psi_soil, ParCost _par_cost, ParPhotosynth _par_photosynth, ParPlant _par_plant, ParEnv _par_env){
    
	const int n = 1;
	// Set up parameters
	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-4;
	param.epsilon_rel = 1e-4;
	param.past = 1;
	param.delta = 5e-5;
	param.max_iterations = 100;

	// Create solver and function object
	LBFGSpp::LBFGSSolver<double> solver(param);
	PHydro_Profit_Inst profit_fun(vcmax, jmax, psi_soil, _par_cost, _par_photosynth, _par_plant, _par_env);

	// bounds
	VectorXd lb(n), ub(n);
	lb << 0.0001;
	ub << 20;

	// Initial guess
	VectorXd x(n);
	x << 1.0; 

	// x will be overwritten to be the best point found
	double fx;
	int niter = solver.minimize(profit_fun, x, fx); //, lb, ub);

	double res = x[0];

	return res;

}


} // phydro

#endif
