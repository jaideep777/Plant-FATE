#ifndef DEMO_RED_MODEL_H
#define DEMO_RED_MODEL_H

#include <individual_base.h>

class LightEnvironment : public EnvironmentBase{
	
	double E = 0;

	public:
	double evalEnv(double x, double t){
		return E;
	}

	// This function must do any necessary precomputations to facilitate evalEnv()
	// Therefore, this should calculate env for all X when it is a function of X
	// In such a case, the solver's SubdivisionSpline can be ussed
	// Note: The state vector in the solver will not be updated until the RK step is completed. 
	// Hence, explicitly pass the state to this function.
	void computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt){
		//             _xm 
		// Calculate _/ w(z,t)u(z,t)dz
		//         xb
		auto w = [S](int i, double t) -> double {
			double z = S->species_vec[0]->getX(i)[0];
			return 0.396*pow(z, 0.749)/10000;
		};
		E = S->state_integral(w, t, 0);
	}

};



class RED_Plant : public IndividualBase <1>{
	public:

	//double input_seed_rain = 1;	
	vector <string> varnames;

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	double a0 = 0.396;
	double phiA = 0.749;
	double g0 = 0.0838;
	double phiG = 0.7134;
	double m0 = 1;
	double mort = 0.035;
	double alpha = 0.1;
	double mu0;

	RED_Plant() {
		mu0 = mort*m0/g0;
	}


	void set_size(const std::array <double, 1>& _x){
	}

	double init_density(void * env, double input_seed_rain){
		return 100/pow(x[0],4);
	}

	void preCompute(double t, void * env){
	}

	double establishmentProbability(double t, void * env){
		return 1;
	}

	std::array<double,1> growthRate(double t, void * env){
		++nrc;
		return {g0*pow(x[0],phiG)};	
	}

	double mortalityRate(double t, void * env){
		++ndc;
		return mort;
	}

	double birthRate(double t, void * env){
		++nbc;
		if (x[0] < 0) throw std::runtime_error("x is negative");
		LightEnvironment* env1 = (LightEnvironment*)env;
		return 0.1/0.9*g0*pow(x[0],phiG)*(1-env1->evalEnv(x[0],t));
	}


};


#endif
